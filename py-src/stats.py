import pysam
import re

def write_eval_file(path, true_single, mismapped, unmapped, multimapped_no_true,
    multimapped_has_true):
    """
    """
    with open(path, "w") as f_out:
        f_out.write("true_single: {}\n".format(true_single))
        f_out.write("mismapped: {}\n".format(mismapped))
        f_out.write("unmapped: {}\n".format(unmapped))
        f_out.write("multimapped_no_true: {}\n".format(multimapped_no_true))
        f_out.write("multimapped_has_true: {}\n".format(multimapped_has_true))
        total = true_single + mismapped + unmapped + multimapped_has_true + \
            multimapped_no_true
        if total != 0:
            f_out.write("mapping rate: {}\n".format( 1 - unmapped / float(total) ) )
            f_out.write("recovered at least 1 correct: {}\n".format( (true_single + multimapped_has_true) / float(total) ) )
        else:
            f_out.write("mapping rate: 0\n")
            f_out.write("recovered at least 1 correct: 0\n" )


def compare_true_observed_positions(true, observed, margin=2):
    """
    """
    return abs(true - observed) < margin

def parse_read_name(name):
    # parts = name.strip("_").split("_")
    # read_id = parts[0]
    # true_position = parts[1]
    pattern = "(\d+)_(\d+)(?:(?:_m=((?:\d+_)+))?)(?:(?:_i=((?:\d+_)+))?)"
    m = re.search(pattern, name)
    if m == None:
        print("Invalid read name -- can not match w/ regex")
        print(name)
        print("Expected format: <id>_<start_pos>(_m=<mismatches>?)(_i=<deletions and ins>?) ")
        print("where mismatches/indesl are a series of postions w/in the read, _ separated")
        exit(1)
    read_id = int(m.group(1))
    true_position = int(m.group(2))
    mismatches = m.group(3)
    indels = m.group(4)
    if mismatches != None:
        mismatches = list(map(int, str(mismatches).strip("_").split("_") ) )
    if indels != None:
        indels = list(map(int, str(indels).strip("_").split("_") ) )

    return true_position, mismatches, indels

def compute_basic_stats(in_path, out_path, margin=3):
    # with open(out_path, "w") as f_out:
    true_single = 0
    mismapped = 0
    f_mismapped = open(out_path + ".mismapped", "w")
    unmapped = 0
    multimapped_has_true = 0
    multimapped_no_true = 0
    with open(in_path, "r") as f_in:
        for line in f_in:
            parts = line.strip().split("\t")

            if len(parts) < 2:
                unmapped += 1
                continue
            true_location, mismatches, indels = parse_read_name(parts[0])

            observed_locations = parts[1]
            observed_locations = list(map(int,observed_locations.split(" ") ) )

            # observed_locations = list(map(int, parts[1:]))

            if len(observed_locations) == 1:
                if compare_true_observed_positions(true_location, observed_locations[0], margin):
                    true_single += 1
                else:
                    f_mismapped.write("{}\n".format(true_location))
                    mismapped += 1
            elif len(observed_locations) == 0:
                unmapped += 1
            else:
                has_true = False
                for obs in observed_locations:
                    if compare_true_observed_positions(obs, true_location, margin):
                        has_true = True
                        break
                if has_true: multimapped_has_true += 1
                else: multimapped_no_true += 1
    write_eval_file(out_path, true_single, mismapped, unmapped, multimapped_no_true, multimapped_has_true)
    # f_out.write("true_single: {}\n".format(true_single))
    # f_out.write("mismapped: {}\n".format(mismapped))
    # f_out.write("unmapped: {}\n".format(unmapped))
    # f_out.write("multimapped_no_true: {}\n".format(multimapped_no_true))
    # f_out.write("multimapped_has_true: {}\n".format(multimapped_has_true))
    # total = true_single + mismapped + unmapped + multimapped_has_true + multimapped_no_true
    # f_out.write("mapping rate: {}\n".format( 1 - unmapped / float(total) ) )
    # f_out.write("recovered at least 1 correct: {}\n".format( (true_single + multimapped_has_true) / float(total) ) )
    f_mismapped.close()

"""
Compute % aligned, aligned uniquely, unmapped, etc for BWA
"""
def compute_basic_bwa_stats(in_path, out_path, margin=3):
    f_mismapped = open(out_path + ".mismapped", "w")
    # parse the BAM file and count stuff
    samfile = pysam.AlignmentFile(in_path, "rb", check_header=False, check_sq=False)
    true_to_observed = {}
    unmapped = 0
    for alignment in samfile.fetch(until_eof=True):
        # true_obs = int(alignment.query_name.split("_")[1])
        read_name = alignment.query_name
        if alignment.is_unmapped:
            unmapped += 1
        else:
            # get the coordinate to which the read aligned
            observed = alignment.reference_start
            if not (read_name in true_to_observed):
                true_to_observed[read_name] = []
            true_to_observed[read_name].append(observed)
    samfile.close()
    # calculate all the metrics
    true_single = 0
    mismapped = 0
    multimapped_has_true = 0
    multimapped_no_true = 0
    margin = 2
    for read_name, values in true_to_observed.items():
        true_position = int(read_name.split("_")[1])
        if len(values) == 1:
            if compare_true_observed_positions(true_position, values[0], margin):
                true_single += 1
            else:
                mismapped += 1
        else:
            # strict comparison
            # if true_position in values:
            found = False
            for v in values:
                if compare_true_observed_positions(true_position, v, margin):
                    found = True
            if found:
                multimapped_has_true += 1
            else:
                multimapped_no_true += 1
    write_eval_file(out_path, true_single, mismapped, unmapped, multimapped_no_true, multimapped_has_true)


if __name__ == "__main__":
    import sys
    import stats
    stats.compute_basic_stats(sys.argv[1], sys.argv[2])
