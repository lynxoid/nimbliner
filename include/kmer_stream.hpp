#ifndef STREAM_KMER_LIB
#define STREAM_KMER_LIB

#include <string>

#include "JellyfishUtil.h"
#include "kseq.h"
#include "definitions.hpp"

namespace nimble {

class KmerStream {
    // TODO: move into the namespace scope
    kmer_t get_all_ones2(const int K) {
        return ( ((kmer_t)1) << ( 2 * K ) ) - (kmer_t)1;
    }

    const char * str;

    size_t str_length;

    bool first;

    kmer_t prev_kmer;

    uint32_t K; // kmer size

    size_t pos = 0;

    kmer_t get_next_kmer(kmer_t bin_kmer, const char next_base, const short K) {
      // mask is all ones
      kmer_t mask = get_all_ones2(K - 1);
      bin_kmer = (bin_kmer & mask) << 2;
      // append a new char on the right
      bin_kmer = bin_kmer | (kmer_t)nimble::dna_codes[ next_base ];
      return bin_kmer;
    }

public:
    KmerStream(const char * s, size_t l, const uint32_t k): str(s), str_length(l), K(k), first(true), pos(0) {
        if (l < K) {
            cerr << "String shorter than the kmer size" << endl;
            exit(1);
        }
    }

    kmer_t getNextBinKmer() {
        if (first) {
            kmer_t bin_kmer = nimble::mer_string_to_binary(str, K);
            kmer_t bin_kmer_next;
            first = false;
            pos++;
            prev_kmer = bin_kmer;
            return bin_kmer;
        }
        else {
            // TODO: if we use 32-mers there may be a collision w/ TTT..TT sequence
            if (pos > str_length - K) return (kmer_t)-1;
            auto next_kmer = get_next_kmer(prev_kmer, str[pos + K - 1], K);
            prev_kmer = next_kmer;
            pos++;
            return next_kmer;
        }
    }
};

}

#endif
