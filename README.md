### Compile

You need libbf. Dont forget about DYLD_LIBRARY_PATH / DYLD_FALLBACK_LIBRARY_PATH.

### Run

To index a chromosome:

```
	./kmers index 20 chromo.fa
```

where 20 is kmer length.

To sample a million reads from the chromosome:

```
	./sample 20 1000000 chromo.fa > sampled_reads.fa
```

To align sampled reads:

```
	./kmers query 20 sampled_reads.fa
```

chro20 genomic reads were aligned to genome assembly at:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

## Summary 

Nimbliner -- rapid read aligner with low memory consumption

References: stored as a pdBG

Reads 1: short, low error rate.

Reads 2: long, high error rate.

## Goals: 
 
- [ ] engineer a testing pipeline; write tests, python code for plotting; script to install/run
  -  [ ] download WGS data for illumina (reference, reads)
  -  [ ] snakemake pipeline
  -  [ ] testing frameworks for python, C++, seq. data
- [ ] integrate with TravisCI
- [ ] find a mapping location for a read, provide mapping quality score, produce positions for reads w/ errors (more kmers sampled)
- [ ] alignment on the laptop (if less memory -- slower, if more memory -- faster)
- [ ] produce full alignments, not just the positions
- [ ] produce full SAM/BAM
- [ ] test speed/performance on long reads

## Comparisons:
 - DALIGNER
 - STAR
 - RapMap
