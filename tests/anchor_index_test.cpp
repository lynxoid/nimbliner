#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <vector>
#include <unordered_map>

#include <bf.h>

// #include "reference_index.hpp"
#include "anchor_index.hpp"
#include "definitions.hpp"
// #include "bloom_reference_index.hpp"

TEST_CASE("AnchorIndex: read anchors", "[readers]")
{
    SECTION("Test hxb2 anchor set")
    {
        nimble::AnchorIndex anchor_idx;
        anchor_idx.readAnchors("data/test/hxb2", 20);
        SECTION("test the number of anchors read")
        {
            REQUIRE(anchor_idx.size() == 194);
        }

        SECTION("test presence/abscence of certain anchors")
        {
            REQUIRE(anchor_idx.is_anchor(152413962224) );
            REQUIRE(!anchor_idx.is_anchor(0) );
        }

        SECTION("test anchor locations")
        {
            vector<seed_position_t> v1 = anchor_idx.get_anchor_locations(152413962224);
            REQUIRE(v1.size() == 1);

            v1 = anchor_idx.get_anchor_locations(912552618219);
            REQUIRE(v1.size() == 2);
            // assert sorted order
            seed_position_t seed = v1[0];
            REQUIRE(seed.first == 0);
            REQUIRE(seed.second == 115);
            seed = v1[1];
            REQUIRE(seed.first == 0);
            REQUIRE(seed.second == 9200);
            v1 = anchor_idx.get_anchor_locations(0);
            REQUIRE(v1.size() == 0);
        }
    }
}
