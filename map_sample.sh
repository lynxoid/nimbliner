#!/bin/bash
set -exou pipefail
# compile w/ docker
IMG=lynxoid/nimbliner_p:0.0.1
# echo "Building fresh docker image"
time docker build -q -t $IMG -f docker/Dockerfile .

SAMPLE="SDBB-NORM-W044216563280-cfDNA-VER-rep1_S12_L001_R1_001"
# time aws s3 cp s3://grail-scna/Verifi/v2/samples/SDBB-NORM-W044216563280-cfDNA-VER-rep1/fastq/HTNVGCCXX/${SAMPLE}.fastq.gz  /mnt/data/bam/${SAMPLE}/

# time gzip -d /mnt/data/bam/${SAMPLE}/${SAMPLE}.fastq.gz

# create genome index
time docker run --rm -t -v /mnt/data/reference/:/reference/ \
			-v /mnt/data/bam/${SAMPLE}/:/data/ \
			$IMG \
			./bin/mapper -i /data/${SAMPLE}.fastq \
                                     -x /reference/decoy_and_viral_nimbliner.0.0.1 \
				     -o /data/${SAMPLE}.nimbliner. \
				     > /mnt/data/bam/${SAMPLE}/${SAMPLE}.nimbliner.alignments
