#include <catch2/catch_test_macros.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <TriangleParametricAtlas.hpp>
#include <TriangleMeshCircleMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>
#include <Tracing.hpp>
#include <Foliation.hpp>
#include <LevelSetCMap.hpp>
#include <ReversedCombinatorialMap.hpp>

double atan2( const Eigen::Vector2d& v )
{
    return atan2( v[1], v[0] );
}

TEST_CASE( "Theta values on cube foliation" )
{
    const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();
    const topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = reparam::sweepEmbedding( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

    const topology::CombinatorialMapBoundary bdry( map );

    const auto bdry_vertex_ids = indexingOrError( bdry, 0 );
    const auto map_face_ids = indexingOrError( map, 2 );

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

    const auto vertex_positions = [&sweep_input]( const topology::CombinatorialMap& map ){
        const auto vertex_ids = indexingOrError( map, 0 );
        return [&sweep_input, vertex_ids]( const topology::Vertex& v ) -> Eigen::Vector3d {
            return sweep_input.mesh.points.at( vertex_ids( v ) );
        };
    };

    const topology::Edge start_edge( phi( bdry, 2, topology::Dart( 0 ) ).value() );
    std::cout << "Right edge? " << keep_face_sides( start_edge.dart() ) << std::endl;
    const reparam::Trace trace = reparam::traceBoundaryField( sides, start_edge, 0.5, sol, vertex_positions( sides ), false );

    const std::vector<double> level_set_values{ 0.0, 0.5, 1.0 };
    const std::vector<reparam::TraceLevelSetIntersection> intersections = reparam::levelSetIntersections( trace, sides, level_set_values );
    REQUIRE( intersections.size() == level_set_values.size() );

    const auto vertex_ids = indexingOrError( map, 0 );

    const auto func = [&]( const topology::Vertex& v ){
        return sol( vertex_ids( v ) );
    };

    { // Base level set
        const auto base_positions = vertex_positions( base );
        const auto face_ids_of_edge = [&]( const topology::Edge& e ){
            return map_face_ids( bdry.toUnderlyingCell( topology::Face( phi( bdry, 2, e.dart() ).value() ) ) );
        };
        const std::map<topology::Vertex, double> thetas = reparam::thetaValues( base, base_positions, face_ids_of_edge, intersections[0] );

        for( const auto& pr : thetas )
        {
            std::cout << "Position: " << base_positions( pr.first ).transpose() << " vs theta: " << pr.second << std::endl;
            const double expected_theta = atan2( base_positions( pr.first ).head<2>() - Eigen::Vector2d( 0.5, 0.5 ) ) + std::numbers::pi / 2;
            CHECK( util::angleEquals( expected_theta, pr.second, 1e-10 ) );
        }
    }

    { // Midway level set
        const topology::LevelSetCMap level_set( map, func, level_set_values[1] );
        const auto v_pos = vertex_positions( map );
        const auto level_set_positions = topology::levelSetVertexPositions( level_set, v_pos );
        const auto face_ids_of_edge = [&]( const topology::Edge& e ){
            return map_face_ids( level_set.underlyingCell( e ) );
        };
        const std::map<topology::Vertex, double> thetas = reparam::thetaValues( level_set, level_set_positions, face_ids_of_edge, intersections[1] );

        for( const auto& pr : thetas )
        {
            std::cout << "Position: " << level_set_positions( pr.first ).transpose() << " vs theta: " << pr.second << std::endl;
            const double expected_theta = atan2( level_set_positions( pr.first ).head<2>() - Eigen::Vector2d( 0.5, 0.5 ) ) + std::numbers::pi / 2;
            CHECK( util::angleEquals( expected_theta, pr.second, 1e-10 ) );
        }
    }

    { // target level set
        const auto target_positions = vertex_positions( target );
        const topology::ReversedCombinatorialMap rev_map( target );
        const auto rev_positions = reversedVertexPositions( rev_map, target_positions );
        const auto face_ids_of_edge = [&]( const topology::Edge& e ){
            return map_face_ids( bdry.toUnderlyingCell( topology::Face( phi( bdry, 2, rev_map.toUnderlyingCell( e ).dart() ).value() ) ) );
        };
        const std::map<topology::Vertex, double> thetas = reparam::thetaValues( rev_map, rev_positions, face_ids_of_edge, intersections[2] );

        for( const auto& pr : thetas )
        {
            std::cout << "Position: " << rev_positions( pr.first ).transpose() << " vs theta: " << pr.second << std::endl;
            const double expected_theta = atan2( rev_positions( pr.first ).head<2>() - Eigen::Vector2d( 0.5, 0.5 ) ) + std::numbers::pi / 2;
            CHECK( util::angleEquals( expected_theta, pr.second, 1e-10 ) );
        }
    }
}
