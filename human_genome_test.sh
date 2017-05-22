set -exou pipefail

# compile w/ docker
TAG=lynxoid/nimbliner:0.2
echo "Building fresh docker image"
docker build -q -t $TAG -f docker/Dockerfile .

DATA="data/hg"
REF="${DATA}/reference"
# create an index
mkdir -p ${DATA}
rm ${REF}.fofn
# echo /nimbliner/${DATA}/chr1.fa > ${REF}.fofn
# echo /nimbliner/${DATA}/chr10.fa >> ${REF}.fofn
echo /nimbliner/${DATA}/chr20.fa >> ${REF}.fofn
CMD="/usr/bin/time -v ./bin/indexer -i ${REF}.fofn -k 15 -o ${REF}"
docker run -i -v $PWD/data/:/nimbliner/data/ $TAG ${CMD}

# align
echo ">read1" > ${DATA}/reads.fa
echo "AACCGTSGTAGTAGTAACA" >> ${DATA}/reads.fa
echo ">read2" >> ${DATA}/reads.fa
echo "TTAAGTAGCATGACGTAGC" >> ${DATA}/reads.fa

CMD="/usr/bin/time -v ./bin/mapper -i ${DATA}/reads.fa -x ${REF} '>' ${REF}.out"
docker run -i -v $PWD/data/:/nimbliner/data/ $TAG ${CMD}
# diff -qiw ${DATA}/reference.expected_output ${REF}.out
