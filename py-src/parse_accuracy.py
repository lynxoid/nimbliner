import sys

no_star = 0
exact_match = 0
one_of_locs = 0
other_locs_no_match = 0;
locs_len = []
with open(sys.argv[1], "r") as f:
	f.readline()
	for line in f:
		locs = map(int, line.strip().split())
		if len(locs) == 1:
			no_star += 1
		if len(locs) == 2:
			if locs[0] == locs[1]:
				exact_match += 1
			else:
				other_locs_no_match += 1
		else:
			locs_len.append( len(locs) - 1 )
			found = False
			for loc in locs[1:]:
				if loc == locs[0]:
					one_of_locs +=1
					found = True
					break
			if not found:
				other_locs_no_match += 1

print "Exact match:", exact_match
print "Not covered by a star:", no_star
print "Was one of the putative locations:", one_of_locs
print "Other locations suggested, original not recovered:", other_locs_no_match
print "Avg number of multimaps (only when >1 locations suggested):", sum(locs_len) * 1.0 / len(locs_len)
			
