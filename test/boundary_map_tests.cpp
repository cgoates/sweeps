#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <SimplexUtilities.hpp>
#include <AbaqusInput.hpp>
#include <Logging.hpp>

using namespace topology;

TEST_CASE( "Importing a 12-tet hex into cmap has correct number of elements", "[single-file]" )
{
    const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();

    const TetMeshCombinatorialMap interior_map( sweep_input.mesh );
    const CombinatorialMapBoundary map( interior_map );

    REQUIRE( cellCount( map, 2 ) == 12 );
    REQUIRE( cellCount( map, 1 ) == 18 );
    REQUIRE( cellCount( map, 0 ) == 8 );

    iterateDartsWhile( map, [&]( const Dart& d ) {
        REQUIRE( not onBoundary( map, d ) );
        return true;
    } );
}

TEST_CASE( "Importing two tets gives correct number of cells", "[single-file]" )
{
    const SweepInput sweep_input = SweepInputTestCases::twoTets();

    const TetMeshCombinatorialMap interior_map( sweep_input.mesh );
    const CombinatorialMapBoundary map( interior_map );

    REQUIRE( cellCount( map, 2 ) == 6 );
    REQUIRE( cellCount( map, 1 ) == 9 );
    REQUIRE( cellCount( map, 0 ) == 5 );

    iterateDartsWhile( map, [&]( const Dart& d ) {
        REQUIRE( not onBoundary( map, d ) );
        return true;
    } );
}

TEST_CASE( "1x1x2 hex each divided into 12 tets", "[single-file]" )
{
    const SweepInput sweep_input = SweepInputTestCases::refinedCube( { 1, 1, 2 } );
    const TetMeshCombinatorialMap interior_map( sweep_input.mesh );
    const CombinatorialMapBoundary map( interior_map );

    REQUIRE( cellCount( map, 2 ) == 20 );
    REQUIRE( cellCount( map, 1 ) == 30 );
    REQUIRE( cellCount( map, 0 ) == 12 );

    iterateDartsWhile( map, [&]( const Dart& d ) {
        REQUIRE( not onBoundary( map, d ) );
        return true;
    } );
}

TEST_CASE( "Simplest gmsh mesh", "[single-file]" )
{
    const SweepInput sweep_input =
        io::loadINPFile( SRC_HOME "/test/simple_mesh.inp", "Surface1", "Surface28" );
    const TetMeshCombinatorialMap interior_map( sweep_input.mesh );
    const CombinatorialMapBoundary map( interior_map );

    REQUIRE( cellCount( map, 2 ) == 24 );
    REQUIRE( cellCount( map, 1 ) == 36 );
    REQUIRE( cellCount( map, 0 ) == 14 );
}