#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>
#include <Logging.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <SimplicialComplexTestCases.hpp>

TEST_CASE( "Laplace patch test", "" )
{
    SweepInput sweep_input = SweepInputTestCases::twelveTetCube();
    sweep_input.mesh.points.back() = Eigen::Vector3d( 0.37, 0.49, 0.55 );

    topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = solveLaplaceSparse( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

    REQUIRE( equals( sol( 8 ), sweep_input.mesh.points.back()( 2 ), sweep_input.mesh.points.back()( 2 ) * 1e-15 ) );
}
