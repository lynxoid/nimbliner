// #define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <vector>
#include <unordered_map>

#include <bf.h>

// #include "anchor_index.hpp"
#include "bloom_reference_index.hpp"

TEST_CASE("BloomReferenceIndex: read index", "[readers]")
{
    SECTION("Test hxb2 anchor set")
    {
        nimble::BloomReferenceIndex idx;
        idx.readIndex("data/test/hxb2");
        SECTION("test that K -- kmer size -- was read correctly")
        {
            REQUIRE(idx.getK() == 20);
        }

        SECTION("test that we read as many kmers as expected")
        {
            REQUIRE(idx.size() == 9084);
        }

    }
}
