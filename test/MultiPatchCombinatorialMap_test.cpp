#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <MultiPatchCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>

using namespace topology;
TEST_CASE( "Simple 2d multipatch cmap" )
{
    const auto cmap_1d_1 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto cmap_1d_2 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto cmap_tp_1 = std::make_shared<const TPCombinatorialMap>( cmap_1d_2, cmap_1d_1 );
    const auto cmap_tp_2 = std::make_shared<const TPCombinatorialMap>( cmap_1d_1, cmap_1d_2 );

    {
        const MultiPatchCombinatorialMap cmap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 4 ) } } } );

        CHECK( cellCount( cmap, 2 ) == 4 );
        CHECK( cellCount( cmap, 1 ) == 12 );
        CHECK( cellCount( cmap, 0 ) == 9 );

        const auto check_phi2 = [&]( const Dart& d, const Dart& expected ) {
            const auto phi2 = phi( cmap, 2, d );
            CHECK( phi2.has_value() );
            if( phi2.has_value() )
            {
                CHECK( phi2.value() == expected );
            }
        };
        check_phi2( 1, 12 );
        check_phi2( 5, 8 );
        check_phi2( 8, 5 );
        check_phi2( 12, 1 );

        CHECK( not phi( cmap, 2, Dart( 0 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 3 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 6 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 7 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 10 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 11 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 13 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 14 ) ).has_value() );
    }

    {
        const MultiPatchCombinatorialMap cmap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 2 ) } } } );

        CHECK( cellCount( cmap, 2 ) == 4 );
        CHECK( cellCount( cmap, 1 ) == 12 );
        CHECK( cellCount( cmap, 0 ) == 9 );

        const auto check_phi2 = [&]( const Dart& d, const Dart& expected ) {
            const auto phi2 = phi( cmap, 2, d );
            CHECK( phi2.has_value() );
            if( phi2.has_value() )
            {
                CHECK( phi2.value() == expected );
            }
        };
        check_phi2( 1, 10 );
        check_phi2( 5, 14 );
        check_phi2( 14, 5 );
        check_phi2( 10, 1 );

        CHECK( not phi( cmap, 2, Dart( 0 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 3 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 6 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 7 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 8 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 11 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 12 ) ).has_value() );
        CHECK( not phi( cmap, 2, Dart( 13 ) ).has_value() );
    }
}

TEST_CASE( "Simple 3d multipatch cmap" )
{
    const auto cmap_1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto cmap_1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto cmap_2d_1 = std::make_shared<const TPCombinatorialMap>( cmap_1d_2, cmap_1d_1 );
    const auto cmap_2d_2 = std::make_shared<const TPCombinatorialMap>( cmap_1d_1, cmap_1d_2 );
    const auto cmap_tp_1 = std::make_shared<const TPCombinatorialMap>( cmap_2d_1, cmap_1d_1 );
    const auto cmap_tp_2 = std::make_shared<const TPCombinatorialMap>( cmap_2d_2, cmap_1d_1 );

    {
        const MultiPatchCombinatorialMap cmap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 19 ) } } } );

        CHECK( cellCount( cmap, 3 ) == 4 );
        CHECK( cellCount( cmap, 2 ) == 20 );
        CHECK( cellCount( cmap, 1 ) == 33 );
        CHECK( cellCount( cmap, 0 ) == 18 );

        const auto check_phi3 = [&]( const Dart& d, const Dart& expected ) {
            const auto phi3 = phi( cmap, 3, d );
            CHECK( phi3.has_value() );
            if( phi3.has_value() )
            {
                CHECK( phi3.value() == expected );
            }
        };
        check_phi3( 1, 67 );
        check_phi3( 2, 70 );
        check_phi3( 3, 69 );
        check_phi3( 4, 68 );
        check_phi3( 25, 91 );
        check_phi3( 26, 94 );
        check_phi3( 27, 93 );
        check_phi3( 28, 92 );

        CHECK( not phi( cmap, 3, Dart( 0 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 20 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 13 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 37 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 23 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 47 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 31 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 24 ) ).has_value() );

        CHECK( not phi( cmap, 3, Dart( 48 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 55 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 51 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 71 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 88 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 79 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 72 ) ).has_value() );
        CHECK( not phi( cmap, 3, Dart( 89 ) ).has_value() );
    }

    {
        const MultiPatchCombinatorialMap cmap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 33 ) } } } );

        CHECK( cellCount( cmap, 3 ) == 4 );
        CHECK( cellCount( cmap, 2 ) == 20 );
        CHECK( cellCount( cmap, 1 ) == 33 );
        CHECK( cellCount( cmap, 0 ) == 18 );
    }
}