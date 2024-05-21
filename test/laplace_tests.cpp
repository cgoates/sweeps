#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>
#include <Logging.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <TriMeshCombinatorialMap.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CommonUtils.hpp>

TEST_CASE( "Laplace patch test", "" )
{
    SweepInput sweep_input = SweepInputTestCases::twelveTetCube();
    sweep_input.mesh.points.back() = Eigen::Vector3d( 0.37, 0.49, 0.55 );

    topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = sweepEmbedding( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

    REQUIRE( util::equals( sol( 8 ), sweep_input.mesh.points.back()( 2 ), sweep_input.mesh.points.back()( 2 ) * 1e-15 ) );
}

TEST_CASE( "Tutte embedding patch test" )
{
    // A perterbation of this:
    // *-------*
    // | \   / |
    // |  \ /  |
    // |   *   |
    // |  / \  |
    // | /   \ |
    // *-------*
    SimplicialComplex mesh{
        { { 0, 1, 4 }, { 0, 4, 2 }, { 4, 1, 3 }, { 4, 3, 2 } },
        { { 0.0, 0.0, 2.0 }, { 1.0, 0.0, 1.2 }, { 0.0, 1.0, 0.7 }, { 1.5, 1.2, 0.2 }, { 0.7, 0.6, 0.0 } } };

    const topology::TriMeshCombinatorialMap map( mesh );

    std::vector<Eigen::Vector2d> constraints{ { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 1.0 } };

    const auto vids = indexingOrError( map, 0 );

    const auto constraints_func = [&]( const topology::Vertex& v ) -> std::optional<Eigen::Vector2d> {
        if( not boundaryAdjacent( map, v ) ) return std::nullopt;
        return constraints.at( vids( v ) );
    };

    const auto vert_positions = [&]( const topology::Vertex& v ){
        return mesh.points.at( vids( v ) );
    };

    const Eigen::MatrixX2d tutte = tutteEmbedding( map, vert_positions, constraints_func, false );

    const Eigen::MatrixX2d expected = ( Eigen::MatrixX2d( 5, 2 ) << 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 0.5 ).finished();


    std::cout << tutte << std::endl;
    REQUIRE( tutte.rows() == 5 );
    for( Eigen::Index i = 0; i < 5; i++ )
    {
        CHECK( tutte( i, 0 ) == expected( i, 0 ) );
        CHECK( tutte( i, 1 ) == expected( i, 1 ) );
    }
}
