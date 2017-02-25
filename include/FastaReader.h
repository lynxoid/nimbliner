#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <zlib.h>
#include <cstdio>
#include "kseq.h"
#include "unistd.h"
#include "fcntl.h"

#ifndef KSEQ
KSEQ_INIT(gzFile, gzread)
//KSEQ_INIT(int, read)
#endif

class FastaReader {
    gzFile fp;
    // int fp; // file handler
    kseq_t *seq;
    int l;

public:

    FastaReader(const char * fname) {
      //fp = fopen(fname, "r"); // TODO: check params
      fp = gzopen(fname, "r");
        if (fp == NULL) {
            printf("Could not open file %s\n", fname);
            return;
        }
        seq = kseq_init(fp);
    }

    bool is_open() { return fp != NULL;}

    ~FastaReader() {
        kseq_destroy(seq);
        //fclose( fp );
        gzclose(fp);
    }

    kseq_t * nextSequence() {
        l = kseq_read(seq);
        if (l < 0) return nullptr;
        else return seq;
    }
};

#endif /* FASTA_READER_H */
