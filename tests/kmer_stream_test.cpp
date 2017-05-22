#include "catch.hpp"

#include <string>

#include "kmer_stream.hpp"
#include "JellyfishUtil.h"
#include "kseq.h"
#include "definitions.hpp"

TEST_CASE("KmerStream: streaming kmers", "[stream_kmers]")
{
    std::string seq("ACGTGCA");
    nimble::KmerStream kmer_stream(seq.c_str(), seq.size(), 3);
    REQUIRE(kmer_stream.getNextBinKmer() == nimble::mer_string_to_binary("ACG", 3) );
    // auto kmer = kmer_stream.getNextBinKmer();
    // cerr << nimble::mer_binary_to_string(kmer, 3);
    REQUIRE(kmer_stream.getNextBinKmer() == nimble::mer_string_to_binary("CGT", 3) );
    REQUIRE(kmer_stream.getNextBinKmer() == nimble::mer_string_to_binary("GTG", 3) );
    REQUIRE(kmer_stream.getNextBinKmer() == nimble::mer_string_to_binary("TGC", 3) );
    REQUIRE(kmer_stream.getNextBinKmer() == nimble::mer_string_to_binary("GCA", 3) );
    REQUIRE(kmer_stream.getNextBinKmer() == (kmer_t)-1 );
}
