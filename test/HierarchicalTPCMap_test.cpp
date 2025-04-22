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

    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 22 ) ) == cmap.dartRanges().toLocalDart( Dart( 22 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 5 ) ) == cmap.dartRanges().toLocalDart( Dart( 5 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 44 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 48 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 52 ) ).second.id() == 4 );
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

    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 18 ) ) == cmap.dartRanges().toLocalDart( Dart( 18 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 6 ) ) == cmap.dartRanges().toLocalDart( Dart( 6 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 96 ) ) == cmap.dartRanges().toLocalDart( Dart( 96 ) ) );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 24 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 112 ) ).second.id() == 4 );
    CHECK( cmap.unrefinedAncestorDartOfCell( Face( 66 ) ).second.id() == 6 );
}

TEST_CASE( "Simplest 3d hierarchical cmap" )
{
    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( 4 );
    const auto tp2d_topo_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_1 );
    const auto tp2d_topo_2 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_2 );
    const auto tp3d_topo_1 = std::make_shared<const TPCombinatorialMap>( tp2d_topo_1, topo1d_2 );
    const auto tp3d_topo_2 = std::make_shared<const TPCombinatorialMap>( tp2d_topo_2, topo1d_4 );

    const HierarchicalTPCombinatorialMap cmap( { tp3d_topo_1, tp3d_topo_2 }, {
        { Volume( 24 ) },
        { Volume( 0 ), Volume( 24 ), Volume( 48 ), Volume( 72 ), Volume( 96 ), Volume( 120 ), Volume( 144 ), Volume( 168 ) }
    } );

    CHECK( cellCount( cmap, 3 ) == 9 );
    CHECK( cellCount( cmap, 2 ) == 41 );
    CHECK( cellCount( cmap, 1 ) == 62 );
    CHECK( cellCount( cmap, 0 ) == 31 );

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

    CHECK( n_darts == 8 * 24 + 40 );
}

TEST_CASE( "3d hierarchical cmap bug" )
{
    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( 4 );
    const auto topo1d_8 = std::make_shared<const CombinatorialMap1d>( 8 );

    const auto topo2d_1 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_1 );
    const auto topo2d_2 = std::make_shared<const TPCombinatorialMap>( topo1d_4, topo1d_2 );
    const auto topo2d_4 = std::make_shared<const TPCombinatorialMap>( topo1d_8, topo1d_4 );

    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo2d_1, topo1d_1 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo2d_2, topo1d_2 );
    const auto tp_topo_4 = std::make_shared<const TPCombinatorialMap>( topo2d_4, topo1d_4 );

    const std::vector<std::vector<topology::Cell>> leaf_elems = [&](){
        const HierarchicalTPCombinatorialMap cmap( { tp_topo_1, tp_topo_2, tp_topo_4 }, { { Volume( 0 ), Volume( 24 ) }, {}, {} } );
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

    const HierarchicalTPCombinatorialMap cmap( { tp_topo_1, tp_topo_2, tp_topo_4 }, leaf_elems );

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

TEST_CASE( "3d hierarchical cmap bug 2" )
{
    const auto three_elem_1d = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto six_elem_1d = std::make_shared<const CombinatorialMap1d>( 6 );
    const auto tp_2d_0 = std::make_shared<const TPCombinatorialMap>( three_elem_1d, three_elem_1d );

    const auto topo2d_1 = std::make_shared<const TPCombinatorialMap>( three_elem_1d, three_elem_1d );
    const auto topo2d_2 = std::make_shared<const TPCombinatorialMap>( six_elem_1d, six_elem_1d );

    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo2d_1, three_elem_1d );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo2d_2, six_elem_1d );

    const std::vector<std::vector<topology::Cell>> leaf_elems = [&](){
        std::vector<topology::Cell> all_first_level;
        iterateCellsWhile( *tp_topo_1, 3, [&]( const topology::Cell& c ) {
            all_first_level.push_back( c );
            return true;
        } );
        const HierarchicalTPCombinatorialMap cmap( { tp_topo_1, tp_topo_2 }, { all_first_level, {} } );
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

    const HierarchicalTPCombinatorialMap cmap( { tp_topo_1, tp_topo_2 }, leaf_elems );

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