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

TIME_CMD = "/usr/bin/time -lp nice "

# one rule to rule them all
rule run_pipeline:
	input:
		expand("{work_dir}/{reference}/k_{K}/plots/main.pdf", 
			work_dir=config["work_dir"], reference=config["reference"],
			K=config["K"])


# Compile a full pdf report
rule compile_pdf_report:
	input:
		# "{work_dir}/{reference}/k_{K}/plots/comparison_table.tex", 
		# expand("{{work_dir}}/{{reference}}/k_{{K}}/plots/{method}_comparison_table.tex", 
			# method=["bwa", "nimbliner"]),
		expand("{{work_dir}}/{{reference}}/k_{{K}}/plots/{method}_comparison_table_all.tex",
			method=["bwa", "nimbliner"]),
		"{work_dir}/{reference}/k_{K}/plots/main.tex"
	output:
		"{work_dir}/{reference}/k_{K}/plots/main.pdf"
	shell:
		"cd {wildcards.work_dir}/{wildcards.reference}/k_{wildcards.K}/plots; pdflatex main; echo {output}"


# copy the clean latex template
rule get_text_template_for_main:
	input:
		"plots/main.tex"
	output:
		"{work_dir}/{reference}/k_{K}/plots/main.tex"
	shell:
		"cp {input} {output}"



# Make a comparison table for reads w/ mismatches
rule make_comparison_table_nimbliner:
	input:
		evals=expand("{{work_dir}}/{{reference}}/k_{{K}}/analysis/nimbliner/sampled_{readlen}_{count}_m={mm_rate}pct_d={d_rate}pct_eval.txt",
					zip,
					readlen=[100 for i in range(7)],
					count=[1000 for i in range(7)],
					mm_rate =[0, 0.5, 2.0, 3.0, 0.0, 0.0, 3.0],
					d_rate	=[0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.2]),
		script="py-src/latex.py"
	output:
		"{work_dir}/{reference}/k_{K}/plots/{method}_comparison_table_all.tex"
	run:
		import latex
		latex.write_table(input.evals, output[0])


# Make a comparison table for reads w/ mismatches
rule make_comparison_table_bwa:
	input:
		evals=expand("{{work_dir}}/{{reference}}/analysis/bwa/sampled_{readlen}_{count}_m={mm_rate}pct_d={d_rate}pct_eval.txt",
					zip,
					readlen=[100 for i in range(7)],
					count=[1000 for i in range(7)],
					mm_rate =[0, 0.5, 2.0, 3.0, 0.0, 0.0, 3.0],
					d_rate	=[0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.2]),
		script="py-src/latex.py"
	output:
		"{work_dir}/{reference}/k_{K}/plots/{method}_comparison_table_all.tex"
	run:
		import latex
		latex.write_table(input.evals, output[0])


# summarize performance results from the output and logs
rule evaluate_nimbliner_alignment:
	input:
		align="{work_dir}/{reference}/k_{K}/alignments/nimbliner/{dataset}.aligned",
		script="py-src/stats.py"
	output:
		"{work_dir}/{reference}/k_{K}/analysis/nimbliner/{dataset}_eval.txt"
	run:
		import stats
		stats.compute_basic_stats(input.align, output[0])



# Align reads to the reference using the index; log the performance
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
		TIME_CMD + "./{input.binary} {wildcards.K} {input.reads} {input.index} {input.anchors} > {output} 2> {log}"



# Run the tool to build an index of the reference sequence
rule build_index:
	input:
		code=["include/reference_index_builder.hpp", "src/sample_reads.cpp"],
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
		TIME_CMD + "./{input.binary} {wildcards.K} {input.ref} 2> {log};"
		# drops all_kmers.txt and star_locations.txt files into the current dir
		"mv all_kmers.txt {output.index};"
		"mv anchors.txt {output.anchors}"



# Compile if any of the input file have changed
rule compile_aligner:
	input:
		"src/kmer_location.cpp",
		"include/aligner.hpp",
		"include/reference_index.hpp"
	output:
		"bin/mapper"
	shell:
		"rm -f {output}; make mapper"


# Compile index builder
rule compile_index_builder:
	input:
		"src/index_builder.cpp",
		"include/reference_index_builder.hpp"
	output:
		"bin/indexer"
	shell:
		"rm -f {output}; make indexer"