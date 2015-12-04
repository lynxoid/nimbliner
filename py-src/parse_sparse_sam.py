import sys

with open(sys.argv[1]) as f:
	i = 1
	for line in f:
		flags, offs, read = line.strip().split()
		if ( int(flags) & 16) == 1:
			# reverse-complement
			read = read[-1:0]
		print ">{}".format(i)
		print read
		i+=1