#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <SimplexUtilities.hpp>
#include <AbaqusInput.hpp>
#include <Logging.hpp>

using namespace topology;

TEST_CASE( "Importing a 12-tet hex into cgogn has correct number of elements", "[single-file]" )
{
    const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();

    const TetMeshCombinatorialMap map( sweep_input.mesh );

    REQUIRE( cellCount( map, 3 ) == 12 );
    REQUIRE( cellCount( map, 2 ) == 30 );
    REQUIRE( cellCount( map, 1 ) == 26 );
    REQUIRE( cellCount( map, 0 ) == 9 );
}

TEST_CASE( "Importing two tets gives correct number of cells", "[single-file]" )
{
    const SweepInput sweep_input = SweepInputTestCases::twoTets();

    const TetMeshCombinatorialMap map( sweep_input.mesh );

    REQUIRE( cellCount( map, 3 ) == 2 );
    REQUIRE( cellCount( map, 2 ) == 7 );
    REQUIRE( cellCount( map, 1 ) == 9 );
    REQUIRE( cellCount( map, 0 ) == 5 );
}

TEST_CASE( "1x1x2 hex each divided into 5 tets", "[single-file]" )
{
    const SweepInput sweep_input = SweepInputTestCases::refinedCube( { 1, 1, 2 } );
    const TetMeshCombinatorialMap map( sweep_input.mesh );

    REQUIRE( cellCount( map, 3 ) == 24 );
    REQUIRE( cellCount( map, 2 ) == 58 );
    REQUIRE( cellCount( map, 1 ) == 47 );
    REQUIRE( cellCount( map, 0 ) == 14 );
}

TEST_CASE( "Simplest gmsh mesh", "[single-file]" )
{
    const SweepInput sweep_input =
        io::loadINPFile( "/Users/caleb/sweeps/attempt-sweep/test/simple_mesh.inp", "Surface1", "Surface28" );
    const TetMeshCombinatorialMap map( sweep_input.mesh );

    REQUIRE( cellCount( map, 3 ) == 24 );
    REQUIRE( cellCount( map, 0 ) == 14 );

    // size_t n_boundary_faces = 0;
    // cgogn::foreach_cell( map, [&]( cgogn::CMap3::Face f ) {
    //     if( cgogn::is_incident_to_boundary( map, f ) ) n_boundary_faces++;
    //     return true;
    // } );
    // REQUIRE( n_boundary_faces == 24 );

    // size_t n_boundary_edge = 0;
    // cgogn::foreach_cell( map, [&]( cgogn::CMap3::Edge e ) {
    //     if( cgogn::is_incident_to_boundary( map, e ) ) n_boundary_edge++;
    //     return true;
    // } );
    // REQUIRE( n_boundary_edge == 36 );
}