#include <catch2/catch_test_macros.hpp>
#include <HierarchicalTPCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

using namespace topology;

TEST_CASE( "Simple hierarchical cmap 1" )
{
    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto topo1d_3 = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto topo1d_6 = std::make_shared<const CombinatorialMap1d>( 6 );
    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_2 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo1d_3, topo1d_6 );

    const HierarchicalTPCombinatorialMap cmap( { tp_topo_1, tp_topo_2 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ), Face( Dart( 12 ) ), Face( 16 ), Face( 20 ), Face( 24 ), Face( 28 ), Face( 32 ) }
    } );

    CHECK( cellCount( cmap, 2 ) == 10 );
    CHECK( cellCount( cmap, 1 ) == 27 );
    CHECK( cellCount( cmap, 0 ) == 18 );

    size_t n_darts = 0;
    iterateDartsWhile( cmap, [&]( const Dart& d ){
        const auto maybe_phi = phi( cmap, {1,-1}, d );
        CHECK( maybe_phi.has_value() );
        if( maybe_phi )
        {
            CHECK( maybe_phi.value() == d );
        }
        n_darts++;
        return true;
    } );

    CHECK( n_darts == 42 );

    CHECK( cmap.unrefinedAncestorDart( Dart( 22 ) ) == cmap.dartRanges().toLocalDart( Dart( 22 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 5 ) ) == cmap.dartRanges().toLocalDart( Dart( 5 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 44 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 48 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 52 ) ).second.id() == 4 );
}

TEST_CASE( "Simple hierarchical cmap 2" )
{
    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( 4 );
    const auto topo1d_8 = std::make_shared<const CombinatorialMap1d>( 8 );
    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_2 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_4 );
    const auto tp_topo_3 = std::make_shared<const TPCombinatorialMap>( topo1d_4, topo1d_8 );

    const HierarchicalTPCombinatorialMap cmap( { tp_topo_1, tp_topo_2, tp_topo_3 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ) },
        { Face( Dart( 40 ) ), Face( Dart( 44 ) ), Face( Dart( 56 ) ), Face( Dart( 60 ) ) }
    } );

    CHECK( cellCount( cmap, 2 ) == 8 );
    CHECK( cellCount( cmap, 1 ) == 23 );
    CHECK( cellCount( cmap, 0 ) == 16 );

    size_t n_darts = 0;
    iterateDartsWhile( cmap, [&]( const Dart& d ){
        const auto maybe_phi = phi( cmap, {1,-1}, d );
        CHECK( maybe_phi.has_value() );
        if( maybe_phi )
        {
            CHECK( maybe_phi.value() == d );
        }
        n_darts++;
        return true;
    } );

    CHECK( n_darts == 36 );

    CHECK( cmap.unrefinedAncestorDart( Dart( 18 ) ) == cmap.dartRanges().toLocalDart( Dart( 18 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 6 ) ) == cmap.dartRanges().toLocalDart( Dart( 6 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 96 ) ) == cmap.dartRanges().toLocalDart( Dart( 96 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 24 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 112 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 66 ) ).second.id() == 6 );
}
