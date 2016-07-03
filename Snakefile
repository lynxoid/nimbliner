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


small_set = {}

big_set = {}

"""
Compile a pdf report
"""
rule compile_pdf_report:
	input:
		"{work_dir}/{reference}/k_{K}/plots/comparison_table.tex", 
		"{work_dir}/{reference}/k_{K}/plots/{method}_comparison_table.tex", 
		"{work_dir}/{reference}/k_{K}/plots/{method}_comparison_table_mismatches.tex",
		"{work_dir}/{reference}/k_{K}/plots/main.tex"
	output:
		expand("{work_dir}/{reference}/k_{K}/plots/main.pdf", 
			work_dir=config["work_dir"], reference=config["reference"],
			K=config["K"])
	shell:
		"cd {wildcards.work_dir}/{wildcards.ref}/k_{wildcards.K}/plots; pdflatex main; echo {output}"


"""
"""
rule get_text_template_for_main:
	input:
		"plots/main.tex"
	output:
		"{work_dir}/{reference}/k_{K}/plots/main.tex"
	shell:
		"cp {input} {output}"


"""
Make a comparison table
"""
rule make_comparison_table:
	input:
		evals=expand("{{work_dir}}/{{reference}}/k_{{K}}/analysis/{{method}}/sampled_{readlen}_{count}_eval.txt",
					zip,
					# work_dir=[config["work_dir"] for i in range(4)],
					# reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 300, 1000],
					count=[100, 1000000, 300000, 10000]
					),
		script="py-src/latex.py"
	output:
		"{work_dir}/{reference}/k_{K}/plots/{method}_comparison_table.tex"
	run:
		import latex
		latex.write_table(input.evals, output[0])


"""
Make a comparison table
"""
rule make_mapping_rate_table:
	input:
		nimble_data=expand("{{work_dir}}/{{reference}}/k_{{K}}/analysis/nimbliner/sampled_{readlen}_{count}_eval.txt",
					zip,
					# work_dir=[config["work_dir"] for i in range(4)],
					# reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 300, 1000],
					count=[100, 1000000, 300000, 10000]
					),
		bwa_data=expand("{{work_dir}}/{{reference}}/k_{{K}}/analysis/bwa/sampled_{readlen}_{count}_eval.txt",
					zip,
					# work_dir=[config["work_dir"] for i in range(4)],
					# reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 300, 1000],
					count=[100, 1000000, 300000, 10000]
					),
		script="py-src/latex.py"
	output:
		"{work_dir}/{reference}/k_{K}/plots/comparison_table.tex"
	run:
		import latex
		latex.write_mapping_rate_table(input.nimble_data, input.bwa_data, output[0])


"""
Make a comparison table for reads w/ mismatches
"""
rule make_comparison_table_mismatches:
	input:
		evals=expand("{{work_dir}}/{{reference}}/k_{{K}}/analysis/{{method}}/sampled_{readlen}_{count}_m={mm_rate}pct_eval.txt",
					zip,
					# work_dir=[config["work_dir"] for i in range(4)],
					# reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 100],
					count=[1000, 1000, 1000],
					mm_rate=[0.5, 1, 2],
					# K=[config["K"] for i in range(4)]
					),
		script="py-src/latex.py"
	output:
		"{work_dir}/{reference}/k_{K}/plots/{method}_comparison_table_mismatches.tex"
	run:
		import latex
		latex.write_table(input.evals, output[0])


"""
"""
rule evaluate_nimbliner_alignment:
	input:
		align="{work_dir}/{reference}/k_{K}/alignments/nimbliner/{dataset}.aligned",
		script="py-src/stats.py"
	output:
		"{work_dir}/{reference}/k_{K}/analysis/nimbliner/{dataset}_eval.txt"
	run:
		import stats
		stats.compute_basic_stats(input.align, output[0])


"""
Align reads to the reference using the index; log the performance
"""
rule align_reads_nimbliner:
	input:
		index="{work_dir}/{reference}/k_{K}/index/{reference}.index",
		anchors="{work_dir}/{reference}/k_{K}/index/{reference}.star",
		reads="{work_dir}/{reference}/sampled/{dataset}.fa",
		binary="bin/mapper"
	output:
		"{work_dir}/{reference}/k_{K}/alignments/nimbliner/{dataset}.aligned"
	log:
		"{work_dir}/{reference}/k_{K}/log/align_{dataset}.log"
	shell:
		"/usr/bin/time -lp ./{input.binary} {wildcards.K} {input.reads} {input.index} {input.anchors} > {output} 2> {log}"


"""
Run the tool to build an index of the reference sequence
"""
rule build_index:
	input:
		code=["include/reference_index.hpp", "src/sample_reads.cpp"],
		ref=expand("{input}/{reference}.fa", 
						input=config["input_dir"], 
						reference=config["reference"]),
		binary="bin/indexer"
	log:
		"{work_dir}/{reference}/log/index_{reference}.log"
	output:
		index="{work_dir}/{reference}/k_{K}/index/{reference}.index",
		anchors="{work_dir}/{reference}/k_{K}/index/{reference}.star"
	shell:
		"./{input.binary} {wildcards.K} {input.ref} 2> {log};"
		# drops all_kmers.txt and star_locations.txt files into the current dir
		"mv all_kmers.txt {output.index};"
		"mv star_locations.txt {output.anchors}"


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

"""
Compile index builder
"""
rule compile_index_builder:
	input:
		"src/index_builder.cpp",
		"include/reference_index.hpp"
	output:
		"bin/indexer"
	shell:
		"rm -f {output}; make indexer"