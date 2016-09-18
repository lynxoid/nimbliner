import sys

diffs = []
with open(sys.argv[1], "r") as f:
	prev = 0
	for l, line in enumerate(f):
		i = int(line.strip())
		if i - prev > 1:
			# print(i-prev)
			diffs.append(i - prev)
		prev = i
		if l % 1000000 == 0:
			print(l, "lines", len(diffs))
print(len(diffs), diffs[:100])