#
# Nimble aligner for short and long reads
# Darya Filippova, 2015-2016
#
# darya.filippova@gmail.com
#

configfile: "experiments.json"

include: "snakes/data.snk"
include: "snakes/bwa_alignments.snk"

import sys
sys.path.append("py-src")


"""
Compile a pdf report
"""
rule compile_pdf_report:
	input:
		expand("{work_dir}/plots/comparison_table.tex", work_dir=config["work_dir"]),
		expand("{work_dir}/plots/{method}_comparison_table.tex", 
			work_dir=config["work_dir"], method=config["methods"]),
		expand("{work_dir}/plots/{method}_comparison_table_mismatches.tex",
			work_dir=config["work_dir"], method=config["methods"]),
		expand("{work_dir}/plots/main.tex", work_dir=config["work_dir"])
	params:
		work_dir=config["work_dir"]
	output:
		expand("{work_dir}/plots/main.pdf", work_dir=config["work_dir"])
	shell:
		"cd {params.work_dir}/plots; pdflatex main; echo {output}"


"""
"""
rule get_text_template_for_main:
	input:
		"plots/main.tex"
	output:
		"{work_dir}/plots/main.tex"
	shell:
		"cp {input} {output}"

# rule get_comparison_tables:
	# input:
		# expand("plots/{method}_comparison_table.tex", method=config["methods"])

"""
Make a comparison table
"""
rule make_comparison_table:
	input:
		evals=expand("{work_dir}/{reference}/analysis/{{method}}/sampled_{readlen}_{count}_eval.txt",
					zip,
					work_dir=[config["work_dir"] for i in range(4)],
					reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 300, 1000],
					count=[100, 1000000, 300000, 10000]
					),
		script="py-src/latex.py"
	output:
		"{work_dir}/plots/{method}_comparison_table.tex"
	run:
		import latex
		latex.write_table(input.evals, output[0])


"""
Make a comparison table
"""
rule make_mapping_rate_table:
	input:
		nimble_data=expand("{work_dir}/{reference}/analysis/nimbliner/sampled_{readlen}_{count}_eval.txt",
					zip,
					work_dir=[config["work_dir"] for i in range(4)],
					reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 300, 1000],
					count=[100, 1000000, 300000, 10000]
					),
		bwa_data=expand("{work_dir}/{reference}/analysis/bwa/sampled_{readlen}_{count}_eval.txt",
					zip,
					work_dir=[config["work_dir"] for i in range(4)],
					reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 300, 1000],
					count=[100, 1000000, 300000, 10000]
					),
		script="py-src/latex.py"
	output:
		"{work_dir}/plots/comparison_table.tex"
	run:
		import latex
		latex.write_mapping_rate_table(input.nimble_data, input.bwa_data, output[0])


"""
Make a comparison table for reads w/ mismatches
"""
rule make_comparison_table_mismatches:
	input:
		evals=expand("{work_dir}/{reference}/analysis/{{method}}/sampled_{readlen}_{count}_m={mm_rate}pct_eval.txt",
					zip,
					work_dir=[config["work_dir"] for i in range(4)],
					reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 100],
					count=[1000, 1000, 1000],
					mm_rate=[0.5, 1, 2]
					),
		script="py-src/latex.py"
	output:
		"{work_dir}/plots/{method}_comparison_table_mismatches.tex"
	run:
		import latex
		latex.write_table(input.evals, output[0])


"""
"""
rule evaluate_nimbliner_alignment:
	input:
		align="{work_dir}/{reference}/alignments/nimbliner/{dataset}.aligned",
		script="py-src/stats.py"
	output:
		"{work_dir}/{reference}/analysis/nimbliner/{dataset}_eval.txt"
	run:
		import stats
		stats.compute_basic_stats(input.align, output[0])


# """
# Align 1 read w/ error to the reference using the index; log the performance
# """
# rule run_read_with1error:
# 	input:
# 		index="{work_dir}/{reference}/index/{reference}.index",
# 		stars="{work_dir}/{reference}/index/{reference}.star",
# 		reads="data/input/test_read_1_error.fa",
# 		binary="bin/mapper"
# 	output:
# 		"{work_dir}/{reference}/alignments/nimbliner/test_read_1_error.aligned"
# 	log:
# 		"{work_dir}/{reference}/log/align_test_read1error.log"
# 	params: K=config["K"]
# 	shell:
# 		"/usr/bin/time -lp ./{input.binary} query {params.K} {input.reads} {input.index} {input.stars} > {output} 2> {log}"


"""
Align reads to the reference using the index; log the performance
"""
rule align_reads:
	input:
		index="{work_dir}/{reference}/index/{reference}.index",
		stars="{work_dir}/{reference}/index/{reference}.star",
		reads="{work_dir}/{reference}/sampled/{dataset}.fa",
		binary="bin/mapper"
	output:
		"{work_dir}/{reference}/alignments/nimbliner/{dataset}.aligned"
	log:
		"{work_dir}/{reference}/log/align_{dataset}.log"
	params: K=config["K"]
	shell:
		"/usr/bin/time -lp ./{input.binary} query {params.K} {input.reads} {input.index} {input.stars} > {output} 2> {log}"


"""
Run the tool to build an index of the reference sequence
"""
rule build_index:
	input:
		code=["include/reference_index.hpp", "src/sample_reads.cpp"],
		ref=expand("{input}/{reference}.fa", 
						input=config["input_dir"], 
						reference=config["reference"]),
		binary="bin/mapper"
	params: K=config["K"]
	log:
		"{work_dir}/{reference}/log/index_{reference}.log"
	output:
		index="{work_dir}/{reference}/index/{reference}.index",
		stars="{work_dir}/{reference}/index/{reference}.star"
	shell:
		"./{input.binary} index {params.K} {input.ref} 2> {log};"
		# drops all_kmers.txt and star_locations.txt files into the current dir
		"mv all_kmers.txt {output.index};"
		"mv star_locations.txt {output.stars}"


"""
Compile if any of the input file have changed
"""
rule compile_aligner:
	input:
		"src/kmer_location.cpp",
		"include/aligner.hpp",
		"include/reference_index.hpp"
	output:
		"bin/mapper"
	shell:
		"rm -f {output}; make mapper"