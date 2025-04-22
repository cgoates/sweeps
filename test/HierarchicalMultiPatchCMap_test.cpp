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

    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 16 ) ) ) == cmap.constituents().at(0)->dartRanges().toLocalDart( Dart( 16 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 6 ) ) ) == cmap.constituents().at(0)->dartRanges().toLocalDart( Dart( 6 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 44 ) ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 48 ) ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 52 ) ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 53 ) ) ).second.id() == 5 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 65 ) ) ).second.id() == 5 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 77 ) ) ).second.id() == 5 );

    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 85 ) ) ).second == cmap.refinementLevels().at( 0 )->toGlobalDart( 1, Dart( 5 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 88 ) ) ).second == cmap.refinementLevels().at( 1 )->toGlobalDart( 1, Dart( 0 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 100 ) ) ).second == cmap.refinementLevels().at( 0 )->toGlobalDart( 1, Dart( 4 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 104 ) ) ).second == cmap.refinementLevels().at( 0 )->toGlobalDart( 1, Dart( 4 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( Dart( 108 ) ) ).second == cmap.refinementLevels().at( 0 )->toGlobalDart( 1, Dart( 4 ) ) );
}

TEST_CASE( "3d multipatch hierarchical cmap bug" )
{
    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( 4 );

    const auto topo2d_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_1 );
    const auto topo2d_2 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_2 );
    const auto topo2d_4 = std::make_shared<const TPCombinatorialMap>( topo1d_4, topo1d_4 );

    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo2d_1, topo1d_1 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo2d_2, topo1d_2 );
    const auto tp_topo_4 = std::make_shared<const TPCombinatorialMap>( topo2d_4, topo1d_4 );

    const auto mp_topo_1 = std::make_shared<const MultiPatchCombinatorialMap>(
        MultiPatchCombinatorialMap( { tp_topo_1, tp_topo_1 }, { { { 0, Dart( 7 ) }, { 1, Dart( 19 ) } } } ) );
    const auto mp_topo_2 = std::make_shared<const MultiPatchCombinatorialMap>(
        MultiPatchCombinatorialMap( { tp_topo_2, tp_topo_2 }, mp_topo_1->connections() ) );
    const auto mp_topo_4 = std::make_shared<const MultiPatchCombinatorialMap>(
        MultiPatchCombinatorialMap( { tp_topo_4, tp_topo_4 }, mp_topo_1->connections() ) );

    const std::vector<std::vector<topology::Cell>> leaf_elems = [&](){
        const HierarchicalMultiPatchCombinatorialMap cmap( { mp_topo_1, mp_topo_2, mp_topo_4 }, { { Volume( 0 ), Volume( 24 ) }, {}, {} } );
        std::vector<std::vector<topology::Cell>> out;
        out.push_back( { Volume( 24 ) } );
        out.push_back( {} );
        out.push_back( {} );
        cmap.iterateChildren( Volume( 0 ), 0, [&]( const Cell& c ) {
            out.at( 1 ).push_back( c );
            return true;
        } );
        const size_t offset = 1;
        cmap.iterateChildren( out.at( 1 ).at( offset ), 1, [&]( const Cell& c ) {
            out.at( 2 ).push_back( c );
            return true;
        } );

        out.at( 1 ).erase( out.at( 1 ).begin() + offset );

        return out;
    }();

    const HierarchicalMultiPatchCombinatorialMap cmap( { mp_topo_1, mp_topo_2, mp_topo_4 }, leaf_elems );

    CHECK( cellCount( cmap, 3 ) == 16 );

    size_t n_darts = 0;
    iterateDartsWhile( cmap, [&]( const Dart& d ){
        for( const auto& phis : std::vector<std::vector<int>>{ {1,-1}, {2, 2}, {3, 3} } )
        {
            const auto maybe_phi = phi( cmap, phis, d );
            if( phis.at( 0 ) != 3 ) CHECK( maybe_phi.has_value() );
            if( maybe_phi )
                CHECK( maybe_phi.value() == d );
        }
        n_darts++;
        return true;
    } );

    CHECK( n_darts == 9 * 24 + 56 + 3 * 40 + 3 * 26 );
}

TEST_CASE( "3d multipatch hierarchical cmap bug 2" )
{
    const auto one_elem_1d = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto two_elem_1d = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto three_elem_1d = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto four_elem_1d = std::make_shared<const CombinatorialMap1d>( 4 );
    const auto six_elem_1d = std::make_shared<const CombinatorialMap1d>( 6 );

    const auto tp_1_1 = std::make_shared<const TPCombinatorialMap>( tensorProductCMapFromComponents({ three_elem_1d, three_elem_1d, one_elem_1d }) );
    const auto tp_1_2 = std::make_shared<const TPCombinatorialMap>( tensorProductCMapFromComponents({ three_elem_1d, three_elem_1d, two_elem_1d }) );

    const auto tp_2_1 = std::make_shared<const TPCombinatorialMap>( tensorProductCMapFromComponents({ six_elem_1d, six_elem_1d, two_elem_1d }) );
    const auto tp_2_2 = std::make_shared<const TPCombinatorialMap>( tensorProductCMapFromComponents({ six_elem_1d, six_elem_1d, four_elem_1d }) );


    const auto mp_1 = std::make_shared<const MultiPatchCombinatorialMap>(
        MultiPatchCombinatorialMap( { tp_1_1, tp_1_2 }, { { { 0, Dart( 5 ) }, { 1, Dart( 0 ) } } } ) );
    const auto mp_2 = std::make_shared<const MultiPatchCombinatorialMap>(
        MultiPatchCombinatorialMap( { tp_2_1, tp_2_2 }, mp_1->connections() ) );

    const std::vector<std::vector<topology::Cell>> leaf_elems = [&](){
        std::vector<topology::Cell> all_first_level;
        iterateCellsWhile( *mp_1, 3, [&]( const topology::Cell& c ) {
            all_first_level.push_back( c );
            return true;
        } );
        const HierarchicalMultiPatchCombinatorialMap cmap( { mp_1, mp_2 }, { all_first_level, {} } );
        std::vector<std::vector<topology::Cell>> out;
        out.push_back( { Volume( 0 ),
                         Volume( 24 ),
                         Volume( 48 ),
                         Volume( 72 ),
                         Volume( 96 ),
                         Volume( 120 ),
                         Volume( 144 ),
                         Volume( 168 ),
                         Volume( 192 ),
                         Volume( 312 ) } );
        const std::vector<topology::Cell> children_of( { Volume( 216 ),
                                                         Volume( 240 ),
                                                         Volume( 264 ),
                                                         Volume( 288 ),
                                                         Volume( 336 ),
                                                         Volume( 360 ),
                                                         Volume( 384 ),
                                                         Volume( 408 ),
                                                         Volume( 432 ),
                                                         Volume( 456 ),
                                                         Volume( 480 ),
                                                         Volume( 504 ),
                                                         Volume( 528 ),
                                                         Volume( 552 ),
                                                         Volume( 576 ),
                                                         Volume( 600 ),
                                                         Volume( 624 ) } );
        out.push_back( {} );
        for( const auto & c : children_of )
        {
            cmap.iterateChildren( c, 0, [&]( const Cell& c ) {
                out.at( 1 ).push_back( c );
                return true;
            } );
        }

        return out;
    }();

    const HierarchicalMultiPatchCombinatorialMap cmap( { mp_1, mp_2 }, leaf_elems );

    iterateDartsWhile( cmap, [&]( const Dart& d ){
        for( const auto& phis : std::vector<std::vector<int>>{ {1,-1}, {2, 2}, {3, 3} } )
        {
            const auto maybe_phi = phi( cmap, phis, d );
            if( phis.at( 0 ) != 3 ) CHECK( maybe_phi.has_value() );
            if( maybe_phi )
                CHECK( maybe_phi.value() == d );
        }
        return true;
    } );
}


