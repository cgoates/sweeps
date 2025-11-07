#include <catch2/catch_test_macros.hpp>
#include <ParametricAtlas.hpp>
#include <MultiPatchParametricAtlas.hpp>
#include <Logging.hpp>

using namespace topology;
using namespace param;

TEST_CASE( "Coordinate Transform 2d" )
{
    // Create a simple multipatch atlas with two copies of a patch glued along different faces
    const auto cmap_1d = std::make_shared<topology::CombinatorialMap1d>( 1 );

    const auto atlas_1d = std::make_shared<param::ParametricAtlas1d>( cmap_1d );

    const auto cmap_2d = std::make_shared<topology::TPCombinatorialMap>( cmap_1d, cmap_1d );
    const auto atlas_2d = std::make_shared<param::TPParametricAtlas>( cmap_2d, atlas_1d, atlas_1d );

    const auto test_with_connection = [&]( const Dart& d1, const Dart& d2, const SmallVector<std::pair<size_t, bool>, 3>& expected_transform ) {

        const auto cmap_mp = std::make_shared<topology::MultiPatchCombinatorialMap>(
            topology::MultiPatchCombinatorialMap( { cmap_2d, cmap_2d },
                                                { { { 0, d1 }, { 1, d2 } } } ) );

        const auto atlas_mp = std::make_shared<param::MultiPatchParametricAtlas>(
            cmap_mp, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_2d, atlas_2d } );
        
        const auto transform = param::coordinateTransform( *atlas_mp, topology::Edge( d1 ) );

        CHECK( transform.size() == expected_transform.size() );
        CHECK( transform.at( 0 ) == expected_transform.at( 0 ) );
        CHECK( transform.at( 1 ) == expected_transform.at( 1 ) );
    };

    test_with_connection( Dart( 0 ), Dart( 0 ), { { 0, false }, { 1, false } } );
    test_with_connection( Dart( 0 ), Dart( 1 ), { { 1, false }, { 0, true } } );
    test_with_connection( Dart( 0 ), Dart( 2 ), { { 0, true }, { 1, true } } );
    test_with_connection( Dart( 0 ), Dart( 3 ), { { 1, true }, { 0, false } } );

    test_with_connection( Dart( 1 ), Dart( 0 ), { { 1, true }, { 0, false } } );
    test_with_connection( Dart( 1 ), Dart( 1 ), { { 0, false }, { 1, false } } );
    test_with_connection( Dart( 1 ), Dart( 2 ), { { 1, false }, { 0, true } } );
    test_with_connection( Dart( 1 ), Dart( 3 ), { { 0, true }, { 1, true } } );

    test_with_connection( Dart( 2 ), Dart( 0 ), { { 0, true }, { 1, true } } );
    test_with_connection( Dart( 2 ), Dart( 1 ), { { 1, true }, { 0, false } } );
    test_with_connection( Dart( 2 ), Dart( 2 ), { { 0, false }, { 1, false } } );
    test_with_connection( Dart( 2 ), Dart( 3 ), { { 1, false }, { 0, true } } );

    test_with_connection( Dart( 3 ), Dart( 0 ), { { 1, false }, { 0, true } } );
    test_with_connection( Dart( 3 ), Dart( 1 ), { { 0, true }, { 1, true } } );
    test_with_connection( Dart( 3 ), Dart( 2 ), { { 1, true }, { 0, false } } );
    test_with_connection( Dart( 3 ), Dart( 3 ), { { 0, false }, { 1, false } } );
}

TEST_CASE( "Coordinate Transform 3d" )
{
    // Create a simple multipatch atlas with two copies of a patch glued along different faces
    const auto cmap_1d = std::make_shared<topology::CombinatorialMap1d>( 1 );

    const auto atlas_1d = std::make_shared<param::ParametricAtlas1d>( cmap_1d );

    const auto cmap_2d = std::make_shared<topology::TPCombinatorialMap>( cmap_1d, cmap_1d );
    const auto atlas_2d = std::make_shared<param::TPParametricAtlas>( cmap_2d, atlas_1d, atlas_1d );

    const auto cmap_3d = std::make_shared<topology::TPCombinatorialMap>( cmap_2d, cmap_1d );
    const auto atlas_3d = std::make_shared<param::TPParametricAtlas>( cmap_3d, atlas_2d, atlas_1d );

    const auto test_with_connection = [&]( const Dart& d1, const Dart& d2, const SmallVector<std::pair<size_t, bool>, 3>& expected_transform ) {

        const auto cmap_mp = std::make_shared<topology::MultiPatchCombinatorialMap>(
            topology::MultiPatchCombinatorialMap( { cmap_3d, cmap_3d },
                                                { { { 0, d1 }, { 1, d2 } } } ) );

        const auto atlas_mp = std::make_shared<param::MultiPatchParametricAtlas>(
            cmap_mp, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_3d, atlas_3d } );
        
        const auto transform = param::coordinateTransform( *atlas_mp, topology::Face( d1 ) );

        CHECK( transform.size() == expected_transform.size() );
        CHECK( transform.at( 0 ) == expected_transform.at( 0 ) );
        CHECK( transform.at( 1 ) == expected_transform.at( 1 ) );
        CHECK( transform.at( 2 ) == expected_transform.at( 2 ) );
    };

    test_with_connection( Dart( 0 ), Dart( 5 ), { { 0, true }, { 1, true }, { 2, true } } );
    test_with_connection( Dart( 0 ), Dart( 6 ), { { 1, false }, { 0, false }, { 2, false } } );
    test_with_connection( Dart( 0 ), Dart( 7 ), { { 1, true }, { 2, true }, { 0, true } } );
    test_with_connection( Dart( 0 ), Dart( 12 ), { {0, true}, {1, false}, {2, false} } );
}