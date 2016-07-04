# write a comparison table in latex! yay
import numpy as np

# parse evaluation script for the number of reads mapped correctly
def parse_eval(path):
	true_single = 0
	multimapT = 0
	multimapF = 0 
	unmapped = 0 
	mismapped = 0
	with open(path, "r") as f:
		for line in f:
			if "true_single" in line:
				true_single = int(line.strip().split(":")[1])
			elif "mismapped" in line:
				mismapped = int(line.strip().split(":")[1])
			elif "unmapped" in line:
				unmapped = int(line.strip().split(":")[1])
			elif "multimapped_no_true" in line:
				multimapF = int(line.strip().split(":")[1])
			elif "multimapped_has_true" in line:
				multimapT = int(line.strip().split(":")[1])
	total = true_single + mismapped + unmapped + multimapT + multimapF
	on_target = true_single + multimapT
	return np.array([on_target, true_single, mismapped, multimapT, multimapF, unmapped, total])

# write a table header w/ the given number of columns and column titles
def write_table_header(f_out, columns):
	c = len(columns)
	f_out.write("\\begin{table}[h!]\n")
	f_out.write("\t\\centering\n")
	f_out.write("\t\\begin{tabular}{l " + " ".join(["r" for i in range(c) ] ) + "}\n")
	f_out.write("\t\\toprule\n")
	# f_out.write("\tDataset & \% on target & trueSingle & mismapped & multimapT & multimapF & unmapped \\\\ \n")
	f_out.write("\t{}\\\\ \n".format( " & ".join(columns) ) )
	f_out.write("\t\\midrule\n")

# write a table footer w/ a given caption
def write_table_footer(f_out, caption):
	f_out.write("\t\\bottomrule\n")
	f_out.write("\t\\end{tabular}\n")
	f_out.write("\t\caption{" + caption + "}\n")
	f_out.write("\end{table}\n")

# write a table comparing mapping performance to latex file
def write_table(inputs, output_path):
	with open(output_path, "w") as f_out:
		write_table_header(f_out, "Dataset & \% on target & trueSingle & falseSingle & trueMulti & falseMulti & unmapped".split(" & "))
		for path in inputs:
			method = path.split("/")[-2]
			name = path.split("/")[-1].rstrip("_eval.txt")
			name = name.replace("_", "\_")
			vals = parse_eval(path)
			values = (vals[:-1] / vals[-1] * 100.0).round(2)
			f_out.write("\t{} & {} \\\\\n".format(name, " & ".join( map(str, values ) ) ) )
		write_table_footer(f_out, "Performance on synthetic data w/o errors for \\textbf{" + method + "}")

# write a table comparing mapping rates across the datasets and aligners
def write_mapping_rate_table(inputs1, inputs2, output):
	f_out = open(output, "w")
	write_table_header(f_out, ["Dataset", "Nimbliner", "BWA"])
	for i in range(len(inputs1)):
		path = inputs1[i]
		name = path.split("/")[-1]
		# make latex happy about the underscores
		name = name.replace("_", "\_")
		vals1 = parse_eval(inputs1[i])
		nimbaliner = round( vals1[0] / vals1[-1] * 100, 2)
		vals2 = parse_eval(inputs2[i])
		bwa = round( vals2[0] / vals2[-1] * 100, 2)
		print(vals2, bwa)
		f_out.write("\t{} & {} & {} \\\\\n".format(name, nimbaliner, bwa ) )
	write_table_footer(f_out, "Another table, totally different")