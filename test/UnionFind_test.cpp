#include <catch2/catch_test_macros.hpp>
#include <UnionFind.hpp>

using namespace util;

constexpr bool aligned = true;

TEST_CASE( "UnionFind basic operations", "[UnionFind]" )
{
    UnionFind uf( 10 );

    REQUIRE( uf.numSets() == 10 );

    uf.unite( 1, 2, aligned );
    REQUIRE( uf.numSets() == 9 );
    REQUIRE( uf.findWithOrientation( 1 ).first == uf.findWithOrientation( 2 ).first );
    REQUIRE( uf.findWithOrientation( 1 ).second == uf.findWithOrientation( 2 ).second );
    uf.unite( 2, 3, not aligned );
    REQUIRE( uf.numSets() == 8 );
    REQUIRE( uf.findWithOrientation( 1 ).first == uf.findWithOrientation( 3 ).first );
    REQUIRE( uf.findWithOrientation( 1 ).second != uf.findWithOrientation( 3 ).second );
    REQUIRE( uf.findWithOrientation( 2 ).first == uf.findWithOrientation( 3 ).first );
    REQUIRE( uf.findWithOrientation( 2 ).second != uf.findWithOrientation( 3 ).second );
}
