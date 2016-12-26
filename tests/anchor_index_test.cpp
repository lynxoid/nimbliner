#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <vector>
#include <unordered_map>

#include <bf.h>

// #include "reference_index.hpp"
#include "anchor_index.hpp"
// #include "bloom_reference_index.hpp"

TEST_CASE("AnchorIndex: read anchors", "[readers]")
{
    SECTION("Test hxb2 anchor set")
    {
        nimble::AnchorIndex anchor_idx;
        anchor_idx.readStarLocations("data/test/hxb2.star", 20);
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
            auto v1 = anchor_idx.get_anchor_locations(152413962224);
            REQUIRE(v1.size() == 1);

            v1 = anchor_idx.get_anchor_locations(912552618219);
            REQUIRE(v1.size() == 2);
            REQUIRE(std::find(v1.begin(), v1.end(), 115) != v1.end() );
            REQUIRE(std::find(v1.begin(), v1.end(), 9200) != v1.end() );

            auto v2 = anchor_idx.get_anchor_locations(0);
            REQUIRE(v2.size() == 0);
        }
    }
}
