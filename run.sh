# compile
rm -f bin/mapper
make nb-mapper
name=data/chr20/results/nimbliner/sampled_100_10000_m=1.0pct_d=0.1pct
alignments=$name.align
summary=$name.summary
debug=$name.debug
time bin/mapper -k 20 -i data/chr20/sampled/sampled_100_10000_m=1.0pct_d=0.1pct.fa -x data/chr20/chr20.flat.index -a data/chr20/anchors.sorted.txt > $alignments 2> $debug
python py-src/stats.py $alignments $summary
