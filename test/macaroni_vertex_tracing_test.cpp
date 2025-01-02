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

TEST_CASE( "Tracing from all the vertices in the macaroni", "[slow]")
{
    const SweepInput sweep_input =
        io::loadINPFile( SRC_HOME "/test/data/macaroni.inp", "Surface3", "Surface4" );

    const topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd ans =
        reparam::sweepEmbedding( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals, reparam::LaplaceEdgeWeights::Cotangent /*FIXME*/ );

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

    const auto keep_face_base = [&]( const topology::Face& f ) {
        return sweep_input.zero_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
                sweep_input.zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
                sweep_input.zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };

    const auto keep_face_target = [&]( const topology::Face& f ) {
        return sweep_input.one_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
                sweep_input.one_bcs.at(
                    bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
                sweep_input.one_bcs.at(
                    bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };

    const topology::CombinatorialMapRestriction sides( bdry, keep_face_sides );
    const topology::CombinatorialMapRestriction base( bdry, keep_face_base );
    const topology::CombinatorialMapRestriction target( bdry, keep_face_target );

    const Eigen::Matrix3Xd grad = gradientsWithBoundaryCorrection( map, sides, ans, normals );

    const auto vertex_positions = [&sweep_input]( const topology::CombinatorialMap& map ){
        const auto vertex_ids = indexingOrError( map, 0 );
        return [&sweep_input, vertex_ids]( const topology::Vertex& v ) -> Eigen::Vector3d {
            return sweep_input.mesh.points.at( vertex_ids( v ) );
        };
    };

    const auto bdry_positions = vertex_positions( bdry );
    topology::GlobalCellMarker traced_vertices( map, 0 );

    iterateCellsWhile( base, 0, [&]( const topology::Vertex& v ) {
        traced_vertices.mark( map, bdry.toUnderlyingCell( v ) );
        return true;
    } );

    // Tracing reverse as a parameterization test
    const auto reverse_ans = -1 * ans;
    const auto reverse_grad = -1 * grad;

    iterateCellsWhile( sides, 0, [&]( const topology::Vertex& v ) {
        if( traced_vertices.isMarked( bdry.toUnderlyingCell( v ) ) ) return true;
        traced_vertices.mark( map, bdry.toUnderlyingCell( v ) );
        CHECK_NOTHROW( reparam::traceBoundaryField( sides, v, 1.0, reverse_ans, bdry_positions, false ) );
        return true;
    } );

    const auto pos = vertex_positions( map );
    iterateCellsWhile( map, 0, [&]( const topology::Vertex& v ) {
        if( traced_vertices.isMarked( v ) ) return true;
        CHECK_NOTHROW( reparam::traceField( map, v, pos( v ), reverse_grad, normals, false ) );
        return true;
    } );
}
