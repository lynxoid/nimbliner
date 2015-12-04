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

