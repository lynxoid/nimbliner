def compute_basic_stats(in_path, out_path):
    with open(out_path, "w") as f_out:
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
                true_location = parts[0]
                observed_locations = parts[1]
                observed_locations = observed_locations.split(" ")
                if len(observed_locations) == 1:
                    if true_location == observed_locations[0]:
                        true_single += 1
                    else:
                        f_mismapped.write("{}\n".format(true_location))
                        mismapped += 1
                elif len(observed_locations) == 0:
                    unmapped += 1
                else:
                    has_true = False
                    for obs in observed_locations:
                        if obs == true_location:
                            has_true = True
                            break
                    if has_true: multimapped_has_true += 1
                    else: multimapped_no_true += 1
        f_out.write("true_single: {}\n".format(true_single))
        f_out.write("mismapped: {}\n".format(mismapped))
        f_out.write("unmapped: {}\n".format(unmapped))
        f_out.write("multimapped_no_true: {}\n".format(multimapped_no_true))
        f_out.write("multimapped_has_true: {}\n".format(multimapped_has_true))
        total = true_single + mismapped + unmapped + multimapped_has_true + multimapped_no_true
        f_out.write("mapping rate: {}\n".format( 1 - unmapped / float(total) ) )
        f_out.write("recovered at least 1 correct: {}\n".format( (true_single + multimapped_has_true) / float(total) ) )
        f_mismapped.close()