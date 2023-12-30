#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <cgogn/core/functions/mesh_info.h>
#include <TempUtils.hpp>
#include <AbaqusInput.hpp>
#include <Logging.hpp>

TEST_CASE( "Importing a 12-tet hex into cgogn has correct number of elements", "[single-file]" ) {
    const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();

    cgogn::CMap3 map;
    mapFromInput( sweep_input, map );

    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Volume>( map ) == 12 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Face>( map ) == 30 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Edge>( map ) == 26 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Vertex>( map ) == 9 );
    REQUIRE( cgogn::is_simplicial( map ) );
    REQUIRE( map.dimension == 3 );
}

TEST_CASE( "Importing two tets gives correct number of cells", "[single-file]" ) {
    const SweepInput sweep_input = SweepInputTestCases::twoTets();

    cgogn::CMap3 map;
    mapFromInput( sweep_input, map );

    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Volume>( map ) == 2 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Face>( map ) == 7 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Edge>( map ) == 9 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Vertex>( map ) == 5 );
    REQUIRE( cgogn::is_simplicial( map ) );
    REQUIRE( map.dimension == 3 );
}

TEST_CASE( "1x1x2 hex each divided into 5 tets", "[single-file]" ) {
    const SweepInput sweep_input = SweepInputTestCases::refinedCube( {1, 1, 2} );
    cgogn::CMap3 map;
    mapFromInput( sweep_input, map );

    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Volume>( map ) == 24 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Face>( map ) == 58 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Edge>( map ) == 47 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Vertex>( map ) == 14 );
    REQUIRE( cgogn::is_simplicial( map ) );
    REQUIRE( map.dimension == 3 );
}

TEST_CASE( "Simplest gmsh mesh", "[single-file]" ) {
    const SweepInput sweep_input = io::loadINPFile( "/Users/caleb/sweeps/attempt-sweep/test/simple_mesh.inp", "Surface1", "Surface28" );
    cgogn::CMap3 map;
    mapFromInput( sweep_input, map );

    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Volume>( map ) == 24 );
    REQUIRE( cgogn::nb_cells<cgogn::CMap3::Vertex>( map ) == 14 );

    size_t n_boundary_faces = 0;
    cgogn::foreach_cell( map, [&]( cgogn::CMap3::Face f ) {
        if( cgogn::is_incident_to_boundary( map, f ) ) n_boundary_faces++;
        return true;
    } );
    REQUIRE( n_boundary_faces == 24 );

    size_t n_boundary_edge = 0;
    cgogn::foreach_cell( map, [&]( cgogn::CMap3::Edge e ) {
        if( cgogn::is_incident_to_boundary( map, e ) ) n_boundary_edge++;
        return true;
    } );
    REQUIRE( n_boundary_edge == 36 );

    REQUIRE( cgogn::is_simplicial( map ) );
    REQUIRE( map.dimension == 3 );
}