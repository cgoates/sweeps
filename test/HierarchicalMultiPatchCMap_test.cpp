#include <catch2/catch_test_macros.hpp>
#include <HierarchicalMultiPatchCombinatorialMap.hpp>
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
    const auto tp_topo_1_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_2 );
    const auto tp_topo_1_2 = std::make_shared<const TPCombinatorialMap>( topo1d_3, topo1d_6 );
    const auto tp_topo_2_1 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_1 );
    const auto tp_topo_2_2 = std::make_shared<const TPCombinatorialMap>( topo1d_6, topo1d_3 );

    const auto mp_topo_1 = std::make_shared<const MultiPatchCombinatorialMap>(
        MultiPatchCombinatorialMap( { tp_topo_1_1, tp_topo_2_1 }, { { { 0, Dart( 1 ) }, { 1, Dart( 4 ) } } } ) );
    const auto mp_topo_2 = std::make_shared<const MultiPatchCombinatorialMap>(
        MultiPatchCombinatorialMap( { tp_topo_1_2, tp_topo_2_2 }, mp_topo_1->connections() ) );

    const HierarchicalMultiPatchCombinatorialMap cmap( { mp_topo_1, mp_topo_2 }, {
        { Face( 4 ), Face( 12 ) },
        { Face( 0 ), Face( 4 ), Face( 8 ), Face( 12 ), Face( 16 ), Face( 20 ), Face( 24 ), Face( 28 ), Face( 32 ),
          Face( 72 ), Face( 76 ), Face( 80 ), Face( 96 ), Face( 100 ), Face( 104 ), Face( 120 ), Face( 124 ), Face( 128 ) }
    } );

    CHECK( cellCount( cmap, 2 ) == 20 );
    CHECK( cellCount( cmap, 1 ) == 52 );
    CHECK( cellCount( cmap, 0 ) == 33 );

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

    CHECK( n_darts == 18 * 4 + 2 * 8 );

    CHECK( cmap.unrefinedAncestorDart( Dart( 16 ) ) == cmap.constituents().at(0)->dartRanges().toLocalDart( Dart( 16 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 6 ) ) == cmap.constituents().at(0)->dartRanges().toLocalDart( Dart( 6 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 44 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 48 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 52 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 53 ) ).second.id() == 5 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 65 ) ).second.id() == 5 );
    CHECK( cmap.unrefinedAncestorDart( Dart( 77 ) ).second.id() == 5 );

    CHECK( cmap.unrefinedAncestorDart( Dart( 85 ) ).second == cmap.refinementLevels().at( 0 )->toGlobalDart( 1, Dart( 5 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 88 ) ).second == cmap.refinementLevels().at( 1 )->toGlobalDart( 1, Dart( 0 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 100 ) ).second == cmap.refinementLevels().at( 0 )->toGlobalDart( 1, Dart( 4 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 104 ) ).second == cmap.refinementLevels().at( 0 )->toGlobalDart( 1, Dart( 4 ) ) );
    CHECK( cmap.unrefinedAncestorDart( Dart( 108 ) ).second == cmap.refinementLevels().at( 0 )->toGlobalDart( 1, Dart( 4 ) ) );
}
