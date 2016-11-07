# compile
rm -f bin/mapper
make nb-mapper
alignments=data/chr20/results/nimbliner/sampled_100_10000_m=3.0pct_d=0.1pct.align
time bin/mapper 20 data/chr20/sampled/sampled_100_10000_m=3.0pct_d=0.1pct.fa data/chr20/chr20.flat.index data/chr20/anchors.sorted.txt > $alignments
python py-src/stats.py $alignments
