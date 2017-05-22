set -exou pipefail

# compile w/ docker
TAG=lynxoid/nimbliner:0.2
echo "Building fresh docker image"
docker build -q -t $TAG -f docker/Dockerfile .

REF="data/smoke/reference"
# # run a smoke test
mkdir -p data/smoke
echo ">seq1" > data/smoke/ref1.fa
echo "AACCGGTT" >> data/smoke/ref1.fa

echo ">seq2" > data/smoke/ref2.fa
echo "TTAAGGCC" >> data/smoke/ref2.fa
echo /nimbliner/data/smoke/ref1.fa > ${REF}.fofn
echo /nimbliner/data/smoke/ref2.fa >> ${REF}.fofn
# build an index for this tiny reference
CMD="./bin/indexer -i ${REF}.fofn -k 3 -o ${REF}"

# align
docker run -i -v $PWD/data/:/nimbliner/data/ $TAG ${CMD}
echo ">read1" > data/smoke/reads.fa
echo "AACC" >> data/smoke/reads.fa
echo ">read2" >> data/smoke/reads.fa
echo "TTAA" >> data/smoke/reads.fa
CMD="./bin/mapper -i data/smoke/reads.fa -x ${REF} > ${REF}.out"
docker run -i -v $PWD/data/:/nimbliner/data/ $TAG ${CMD}
diff -qiw data/smoke/reference.expected_output ${REF}.out
