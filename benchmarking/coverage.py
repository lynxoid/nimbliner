# parse anchors and compute coverage for chromo 20
from Bio import SeqIO
import numpy as np

anchors = "data/chr20/index/infreq_anchors/chr20.infreq_anchors.txt"

for record in SeqIO.parse("/Users/lynxoid/dev/seq_data/hs_alt_CHM1_1.1_chr20.fa", "fasta"):
    N = len(record.seq)
    print("Seq length: {}bp".format(N))

N = 1000
# covered locations
L = []
with open(anchors, "r") as f:
    i = 0
    for line in f:
        parts = line.strip().split()
        locations = map(int, parts[1:])
        L += locations
        # print(L)
        i += 1
        if i >= 10: break

print("Covered locations:", len(L))
L = sorted(L)
# print(L[1000:1010])
cov = {i: 0 for i in range(N)}
for l in L:
    for i in range(40):
        cov[l + i - 19] = 1
    # print("{}: [{}, {}]".format(l, l - 19, l+i - 19) )
print("Updated coverage dict")

covered_bases = sum(cov.values())
print("Covered bases", covered_bases)

print("Computing uncovered stretches")
uncovered = []
# compute uncovered stretches
# 0 1 1 0 0 1
# 1 0 0 0 0 0
start = -1
end = -1

mn = min(cov.keys())
mx = min(cov.keys())
for i in range(mn - 1, N):
# for i in range(100):
    if cov[i] == 0:
        if start >= 0:
            end = i
        else:
            start = i
            end = i
        # print("branch 0", i, start, end)
    else:
        print("branch 1", i, cov[i], start, end)
        if start > 0:
            uncovered.append( end - start + 1 )
            start = -1
            end = -1

uncovered.append(end - start + 1)
print(uncovered)
counts, bins = np.histogram(uncovered)
print(counts, bins)
