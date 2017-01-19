## nimbliner - nimble aligner

Fast and lightweight read aligner (experimental) ![CircleCI](https://circleci.com/gh/lynxoid/nimbliner/tree/dev.svg?style=shield&circle-token=:circle-token)

###

Nimbliner uses Bloom filters instead of suffix arrays as reference which incurs the cost close to `n` in the size of reference sequence (instead of `2-4n` for suffix arrays or BWT). It also does not need to perform a full alignment shaving off a lot of the computational cost. Nimbliner does not yet produce cigar strings, but there is no reason why it would not be able to.

### Installation

You can build a docker image and then run nimbliner wthin the container:

```
<clone the repo>
docker build -t nimbliner-dev:0.1 -f docker/Dockerfile docker/
# this will produce all_kmers.txt and anchors.txt in the current directory
docker run -v `pwd`:/nimbliner nimbliner-dev:0.1 indexer 20 <path to your reference, single fasta file>
docker run -v `pwd`:/nimbliner nimbliner-dev:0.1 mapper 20 <path to your reads, single fasta file> all_kmers.txt anchors.txt
```

#### Other ways

You can compile from source. The dependencies are [liffb](https://github.com/mavam/libbf) and [TCLAP](http://tclap.sourceforge.net/). You may need to set `LD_LIBRARY_PATH` (or `DYLD_LIBRARY_PATH` for MacOS) to `/usr/local/lib` since `libbf` installs there by default.

### Niceties

You can generate synthetic reads w/ mismatches and indels. For example, to sample a million reads from the chromosome w/ 1.5% error rate, do:

```
docker run -v `pwd`:/nimbliner nimbliner-dev:0.1 sample 20 1000000 chromo.fa 1.5 > sampled_reads.fa
```

### TODO
- [ ] provide benchmark data
- [ ] prepare indices for the whole human genome
- [ ] integrate with TravisCI

### Comparisons:
 - DALIGNER
 - STAR
 - BWA
 - RapMap (speed-only, RapMap does not generate full alignments)
