#include <catch2/catch_test_macros.hpp>
#include <Tracing.hpp>
#include <Logging.hpp>
#include <SimplexUtilities.hpp>
#include <CombinatorialMap.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>
#include <Laplace.hpp>
#include <AbaqusInput.hpp>
#include <random>

bool trueOnARandomQuarterOfCalls()
{
    static auto gen = std::bind(std::uniform_int_distribution<>(0,4),std::default_random_engine());
    return gen() < 1;
}

TEST_CASE( "Tracing from all the interior faces in the macaroni", "[slow]")
{
    const SweepInput sweep_input =
        io::loadINPFile( SRC_HOME "/test/data/macaroni.inp", "Surface3", "Surface4" );

    const topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd ans =
        solveLaplaceSparse( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

    const topology::CombinatorialMapBoundary bdry( map );

    const auto bdry_vertex_ids = indexingOrError( bdry, 0 );

    const auto keep_face_sides = [&]( const topology::Face& f ) {
        return ( not sweep_input.zero_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) or
                    not sweep_input.zero_bcs.at(
                        bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) or
                    not sweep_input.zero_bcs.at(
                        bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) ) ) and
                ( not sweep_input.one_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) or
                    not sweep_input.one_bcs.at(
                        bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) or
                    not sweep_input.one_bcs.at(
                        bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) ) );
    };

    const topology::CombinatorialMapRestriction sides( bdry, keep_face_sides );

    const Eigen::Matrix3Xd grad = gradientsWithBoundaryCorrection( map, sides, ans, normals );

    size_t i = 0;
    iterateCellsWhile( map, 2, [&]( const topology::Face& start_face ) {
        if( onBoundary( map, start_face.dart() ) ) return true;
        if( i++ % 4 != 0 ) return true; // To shorten the test
        CHECK_NOTHROW( traceField( map, start_face, centroid( map, start_face ), grad, normals ) );
        return true;
    } );
}
