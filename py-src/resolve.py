# sort stuff

A = [2, 18, 25]
B = [7, 19, 22, 30]
x = 5
k = 3

i = 0
j = 0
l_a = len(A)
l_b = len(B)
while (i < l_a and j < l_b):
	if A[i] < B[j]:
		if B[j] - A[i] == x:
			print A[i]
			i += 1
			j += 1
		elif B[j] - A[i] < x:
			# too close -- try moving B[j]
			j += 1
		else:
			# too far apart -- try moving A[i]
			i += 1
	else:
		# A[i] >= B[j]
		j += 1