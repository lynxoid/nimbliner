# # compile
# rm -f bin/mapper
# make nb-mapper
# # echo $PWD
# # ls .
# dir=data/chr20
# name=$dir/results/nimbliner/sampled_100_10000_m=1.0pct_d=0.1pct
# ls -l $name.*
# alignments=$name.align
# summary=$name.summary
# debug=$name.debug
# # bin/sample -d 0.1 -m 1 -n 10000 -l 100 -i /data/human/GRCh38/chr20.fna -o data/chr20/sampled/sampled_100_10000_m=1.0pct_d=0.1pct.fa
# time bin/mapper -k 20 -i $dir/sampled/sampled_100_10000_m=1.0pct_d=0.1pct.fa -x $dir/chr20.flat.index -a $dir/chr20.anchors.txt > $alignments 2> $debug
# python py-src/stats.py $alignments $summary

# compile w/ docker
TAG=lynxoid/nimbliner:0.2
docker build -t $TAG -f docker/Dockerfile .
# run tests
docker run -it -v $PWD/data/:/nimbliner/data/ $TAG;./nimbliner/test
