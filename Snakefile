#
# Nimble aligner for short and long reads
# Darya Filippova, 2015-2016
#
# darya.filippova@gmail.com
#

# config: "experiments.json"

import sys
sys.path.append("py-src")


"""
Compile a pdf report
"""
rule compile_pdf_report:
	input:
		"plots/comparison_table.tex",
		"plots/main.tex"
	output:
		"plots/main.pdf"
	shell:
		"cd plots; pdflatex main"


rule get_comparison_tables:
	input:
		expand("plots/{method}_comparison_table.tex", method=config["methods"])

"""
Make a comparison table
"""
rule make_comparison_table:
	input:
		evals=expand("{work_dir}/{reference}/analysis/{method}}/sampled_{readlen}_{count}_eval.txt",
					zip,
					work_dir=[config["work_dir"] for i in range(4)],
					reference=[config["reference"] for i in range(4)],
					readlen=[100, 100, 300, 1000],
					count=[100, 1000000, 300000, 10000]
					),
		script="py-src/latex.py"
	output:
		"plots/{method}_comparison_table.tex"
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
		"plots/comparison_table.tex"
	run:
		import latex
		latex.write_mapping_rate_table(input.nimble_data, input.bwa_data, output[0])


"""
"""
rule evaluate_alignment:
	input:
		align="{work_dir}/{reference}/alignments/nimbliner/{dataset}.aligned",
		script="py-src/stats.py"
	output:
		"{work_dir}/{reference}/analysis/nimbliner/{dataset}_eval.txt"
	run:
		import stats
		stats.compute_basic_stats(input.align, output[0])


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
Run the tool to build an index
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
"""
rule evaluate_bwa_alignments:
	input:
		align="{work_dir}/{reference}/alignments/bwa/{dataset}.bam",
		script="py-src/stats.py"
	output:
		"{work_dir}/{reference}/analysis/bwa/{dataset}_eval.txt"
	run:
		import stats
		stats.compute_basic_bwa_stats(input.align, output[0])


"""
Align reads using bwa mem
"""
rule align_reads_bwa:
	input:
		index="{work_dir}/{reference}/index/bwa/{reference}.ann",
		reads="{work_dir}/{reference}/sampled/{dataset}.fa"
	params:
		index="{work_dir}/{reference}/index/bwa/{reference}",
	output:
		"{work_dir}/{reference}/alignments/bwa/{dataset}.bam"
	log:
		"{work_dir}/{reference}/log/bwa/align_{dataset}.log"
	params: K=config["K"]
	shell:
		"/usr/bin/time -lp bwa mem {params.index} {input.reads} 2> {log} | samtools view -Sb - > {output} "

"""
Build index for bwa
"""
rule build_bwa_index:
	input:
		ref=expand("{input}/{reference}.fa", 
						input=config["input_dir"], 
						reference=config["reference"]),
	params:
		dir="{work_dir}/{reference}/index/bwa/{reference}",
	output:
		"{work_dir}/{reference}/index/bwa/{reference}.amb",
		"{work_dir}/{reference}/index/bwa/{reference}.ann",
		"{work_dir}/{reference}/index/bwa/{reference}.bwt",
		"{work_dir}/{reference}/index/bwa/{reference}.pac",
		"{work_dir}/{reference}/index/bwa/{reference}.sa"
	shell:
		"bwa index -p {params.dir} {input.ref}"


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
Run the sampling w/o errors
"""
rule sample_reference_no_error:
	input:
		ref=expand("{input}/{reference}.fa", 
					input=config["input_dir"], 
					reference=config["reference"]),
		binary="bin/sample"
	output:
		"{work_dir}/{reference}/sampled/sampled_{R}_{N}.fa"
	log:
		"{work_dir}/{reference}/log/sample_{R}_{N}.log"
	shell:
		"./{input.binary} {wildcards.R} {wildcards.N} {input.ref} > {output} 2> {log};"
		# "ln -s {output} {wilcards.input}/sampled_{wildcards.R}_{wildcards.N}.fa"




"""
Compile if any of the input file have changed
"""
rule compile_sampler:
	input:
		"src/sample_reads.cpp",
	output:
		"bin/sample"
	shell:
		"rm -f {output}; make sample"


"""
Download PacBio lambda virus data
"""
rule download_pacbio_data:
	output:
	shell:
		"echo 'Download pacbio data'"

"""
Download illumina NA128... WGS dataset
"""
rule download_illumina_wgs_data:
	output:
		expand("")
	shell:
		"touch {output}"

"""
Download human chromosome (parametrized)
"""
rule download_human_chr:
	params: "ftp:long_url."
	shell:
		"echo 'download human chr'"