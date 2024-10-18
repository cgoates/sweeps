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
#include <DelaunayTriangulation.hpp>
#include <VTKOutput.hpp>
#include <AbaqusInput.hpp>
#include <iomanip>
#include <memory>

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
            const double expected_theta = atan2( rev_positions( pr.first ).head<2>() - Eigen::Vector2d( 0.5, 0.5 ) ) + std::numbers::pi / 2;
            CHECK( util::angleEquals( expected_theta, pr.second, 1e-10 ) );
        }
    }
}

void testLevelSetBasedTracing( const SweepInput& sweep_input,
                               const std::vector<double> level_set_values,
                               const std::optional<std::string> output_prefix )
{
    reparam::levelSetBasedTracing( sweep_input, level_set_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
        SimplicialComplex level_sets;
        SimplicialComplex param_out;

        if( output_prefix )
        {
            for( const auto& leaf : leaves )
            {
                iterateCellsWhile( leaf.space_mapping->parametricAtlas().cmap(), 2, [&]( const topology::Face& f ) {
                    addTriangleNoDuplicateChecking( level_sets,
                                                    triangleOfFace<3>( leaf.space_mapping->parametricAtlas().cmap(),
                                                                       leaf.space_mapping->vertPositions(),
                                                                       f ) );
                    return true;
                } );
            }
        }

        // Draw some lines!
        param_out.points.reserve( 8 * 6 * level_set_values.size() );
        param_out.simplices.reserve( 8 * 6 * level_set_values.size() );
        for( const double radius : { 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 } )
        {
            for( const double theta : { 0, 45, 90, 135, 180, 225, 270, 315 } )
            {
                const Eigen::Vector2d circle_pt( radius * cos( theta * std::numbers::pi / 180 ),
                                                radius * sin( theta * std::numbers::pi / 180 ) );
                const size_t offset = param_out.points.size();
                size_t i = 0;
                for( const auto& leaf : leaves )
                {
                    const auto& param_pt = leaf.circle_mapping->maybeInverse( circle_pt );
                    CHECK( param_pt.has_value() );
                    if( not param_pt.has_value() )
                    {
                        std::cerr << "NO VALUE r: " << radius << " t: " << theta << " i: " << i << std::endl;
                        break;
                    }
                    if( output_prefix )
                    {
                        const auto space_pt = leaf.space_mapping->evaluate( param_pt.value().first, param_pt.value().second );
                        param_out.points.push_back( space_pt );
                    }
                    i++;
                }

                if( output_prefix )
                {
                    if( i == 0 ) continue;
                    for( size_t simplex_ii = 0; simplex_ii < i - 1; simplex_ii++ )
                    {
                        param_out.simplices.push_back( Simplex( offset + simplex_ii, offset + simplex_ii + 1 ) );
                    }
                }
            }
        }

        if( output_prefix )
        {
            io::VTKOutputObject output( level_sets );
            io::outputSimplicialFieldToVTK( output, output_prefix.value() + "_level_sets.vtu" );

            io::VTKOutputObject output2( param_out );
            io::outputSimplicialFieldToVTK( output2, output_prefix.value() + "_foliation_param.vtu" );
        }
    } );
}

void testLevelSetBasedTracing( const SweepInput& sweep_input,
                               const size_t n_levels,
                               const std::optional<std::string> output_prefix = std::nullopt )
{
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    testLevelSetBasedTracing( sweep_input, level_set_values, output_prefix );
}

TEST_CASE( "Level set parameterization of left ventricle" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/left_ventricle.inp", "Surface3", "Surface2" );

    using util::linspace;
    using util::concatenate;

    // level set every 0.02 until 0.75, then ~30 levels between 0.75 and 0.85, then every 0.2 until 1.0.
    const std::vector<double> level_set_values =
        concatenate(
            concatenate( linspace( 0, 0.78, 40 ), linspace( 0.8, 0.819, 20 ) ),
            concatenate( concatenate( linspace( 0.81925, 0.81975, 3 ), linspace( 0.82, 0.85, 31 ) ), linspace( 0.86, 1.0, 13 ) ) );

    const std::optional<std::string> output_prefix = std::nullopt;//{ "ventricle" };
    testLevelSetBasedTracing( sweep_input, level_set_values, output_prefix );
}

TEST_CASE( "Level set parameterization of hook" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/hook.inp", "Surface12", "Surface10" );

    const std::optional<std::string> output_prefix = std::nullopt;//{ "hook" };
    testLevelSetBasedTracing( sweep_input, 30, output_prefix );
}

TEST_CASE( "Level set parameterization of cube" )
{
    const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();

    testLevelSetBasedTracing( sweep_input, 3 );
}