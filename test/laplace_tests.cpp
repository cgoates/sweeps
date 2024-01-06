#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>

TEST_CASE( "Laplace patch test", "[single-file]" ) {
    SweepInput sweep_input = SweepInputTestCases::twelveTetCube();
    sweep_input.points.back() = Eigen::Vector3d( 0.37, 0.49, 0.55 );

    cgogn::CMap3 map;
    SimplexUtilities::mapFromInput( sweep_input, map );    
    const Eigen::VectorXd sol = solveLaplaceSparse( map, sweep_input.zero_bcs, sweep_input.one_bcs );

    // NOTE: I don't love that this is only within five percent, but I'm not sure if I should expect better...
    // At least this is a smoke test for now.
    REQUIRE( std::abs( sol( 8 ) - sweep_input.points.back()( 2 ) ) < sweep_input.points.back()( 2 ) * 0.05 );
}

