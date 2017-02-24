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
echo "Building fresh docker image"
docker build -q -t $TAG -f docker/Dockerfile .
# run tests
# docker run -i -v $PWD/data/:/nimbliner/data/ $TAG ./bin/nb_tests

# # run a smoke test
mkdir -p data/smoke
echo ">seq1" > data/smoke/ref1.fa
echo "AACCGGTT" >> data/smoke/ref1.fa

echo ">seq2" > data/smoke/ref2.fa
echo "TTAAGGCC" >> data/smoke/ref2.fa
# echo $PWD/data/smoke/ref1.fa > data/smoke/reference.fofn
# echo $PWD/data/smoke/ref2.fa >> data/smoke/reference.fofn
echo /nimbliner/data/smoke/ref1.fa > data/smoke/reference.fofn
echo /nimbliner/data/smoke/ref2.fa >> data/smoke/reference.fofn

docker run -i -v $PWD/data/:/nimbliner/data/ $TAG ./bin/indexer -i data/smoke/reference.fofn -k 3 -o data/smoke/reference
echo ">read1" > data/smoke/reads.fa
echo "AACC" >> data/smoke/reads.fa
echo ">read2" >> data/smoke/reads.fa
echo "TTAA" >> data/smoke/reads.fa
docker run -i -v $PWD/data/:/nimbliner/data/ $TAG ./bin/mapper -i data/smoke/reads.fa -x data/smoke/reference > data/smoke/reference.out
diff -qiw data/smoke/reference.expected_output data/smoke/reference.out
