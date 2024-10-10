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

void testLevelSetBasedTracing( const SweepInput& sweep_input, const std::vector<double> level_set_values, const bool log, const std::optional<std::string> output_prefix )
{
    const topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = reparam::sweepEmbedding( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

    if( log ) std::cout << "FINISHED LAPLACE\n\n";

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
               sweep_input.one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
               sweep_input.one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };

    const topology::CombinatorialMapRestriction sides( bdry, keep_face_sides );
    const topology::CombinatorialMapRestriction base( bdry, keep_face_base, true );
    const topology::CombinatorialMapRestriction target( bdry, keep_face_target, true );

    const auto vertex_positions = [&sweep_input]( const topology::CombinatorialMap& map ) {
        const auto vertex_ids = indexingOrError( map, 0 );
        return [&sweep_input, vertex_ids]( const topology::Vertex& v ) -> Eigen::Vector3d {
            return sweep_input.mesh.points.at( vertex_ids( v ) );
        };
    };

    const topology::Edge start_edge = [&]() {
        std::optional<topology::Edge> out;
        iterateDartsWhile( sides, [&]( const topology::Dart& d ) {
            const auto phi2 = phi( bdry, 2, d );
            if( keep_face_base( phi2.value() ) )
            {
                out.emplace( d );
                return false;
            }
            return true;
        } );
        if( not out.has_value() ) throw std::runtime_error( "Unable to find tracing start edge" );
        return out.value();
    }();
    const reparam::Trace trace =
        reparam::traceBoundaryField( sides, start_edge, 0.5, sol, vertex_positions( sides ), false );

    if( log ) std::cout << "FINISHED TRACE\n\n";

    const std::vector<reparam::TraceLevelSetIntersection> intersections =
        reparam::levelSetIntersections( trace, sides, level_set_values );
    REQUIRE( intersections.size() == level_set_values.size() );

    const auto vertex_ids = indexingOrError( map, 0 );

    const auto harmonic_func = [&]( const topology::Vertex& v ) { return sol( vertex_ids( v ) ); };

    SimplicialComplex level_sets;

    SimplicialComplex param_out;

    struct FoliationLeaf
    {
        // For the middle leaves
        std::unique_ptr<topology::LevelSetCMap> level_set_cmap;
        std::unique_ptr<topology::DelaunayTriangulation> level_set_tri;

        // For the target surface
        std::unique_ptr<topology::ReversedCombinatorialMap> reversed_cmap;

        // For all leaves
        std::unique_ptr<const Eigen::MatrixX2d> tutte;
        std::unique_ptr<param::TriangleParametricAtlas> atlas;
        std::unique_ptr<mapping::TriangleMeshCircleMapping> circle_mapping;
        std::unique_ptr<mapping::TriangleMeshMapping> space_mapping;
    };

    const auto process_param =
        [&]( const topology::CombinatorialMap& cmap, const auto& positions, const auto& thetas, FoliationLeaf& leaf ) {
            const auto base_vert_ids = indexingOrError( cmap, 0 );
            const std::map<size_t, double> thetas_by_id = [&]() {
                std::map<size_t, double> out;
                for( const auto& pr : thetas )
                {
                    out.insert( { base_vert_ids( pr.first ), pr.second } );
                }
                return out;
            }();

            const auto constraints_func = [&]( const topology::Vertex& v ) -> std::optional<Eigen::Vector2d> {
                if( boundaryAdjacent( cmap, v ) )
                {
                    const double theta = thetas_by_id.at( base_vert_ids( v ) );
                    return Eigen::Vector2d( cos( theta ), sin( theta ) );
                }
                else
                    return std::nullopt;
            };

            leaf.tutte = std::make_unique<const Eigen::MatrixX2d>(
                reparam::tutteEmbedding( cmap, positions, constraints_func, true ) );

            leaf.atlas = std::make_unique<param::TriangleParametricAtlas>( cmap );
            const auto vert_positions = [tutte = *( leaf.tutte ), base_vert_ids]( const topology::Vertex& v ) {
                return tutte.row( base_vert_ids( v ) );
            };
            leaf.circle_mapping = std::make_unique<mapping::TriangleMeshCircleMapping>( *leaf.atlas, vert_positions );
            leaf.space_mapping = std::make_unique<mapping::TriangleMeshMapping>( *leaf.atlas, positions, 3 );

            if( output_prefix )
                iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
                    addTriangleNoDuplicateChecking( level_sets, triangleOfFace<3>( cmap, positions, f ) );
                    return true;
                } );
        };

    std::vector<FoliationLeaf> leaves;
    leaves.reserve( level_set_values.size() );

    { // Base level set
        const auto base_positions = vertex_positions( bdry );
        const auto face_ids_of_edge = [&]( const topology::Edge& e ) {
            return map_face_ids( bdry.toUnderlyingCell( topology::Face( phi( bdry, 2, e.dart() ).value() ) ) );
        };
        const std::map<topology::Vertex, double> thetas =
            reparam::thetaValues( base, base_positions, face_ids_of_edge, intersections[0] );

        leaves.push_back( {} );
        process_param( base, base_positions, thetas, leaves.back() );
    }

    if( log ) std::cout << "FINISHED BASE\n\n";

    for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
    { // Midway level set
        leaves.push_back( {} );
        leaves.back().level_set_cmap =
            std::make_unique<topology::LevelSetCMap>( map, harmonic_func, level_set_values[level_ii] );
        const auto v_pos = vertex_positions( map );
        const auto& level_set = *leaves.back().level_set_cmap;
        const auto level_set_positions = topology::levelSetVertexPositions( level_set, v_pos );
        const auto face_ids_of_edge = [&]( const topology::Edge& e ) {
            return map_face_ids( level_set.underlyingCell( e ) );
        };
        const std::map<topology::Vertex, double> thetas =
            reparam::thetaValues( level_set, level_set_positions, face_ids_of_edge, intersections[level_ii] );

        leaves.back().level_set_tri =
            std::make_unique<topology::DelaunayTriangulation>( level_set, level_set_positions );
        const auto tri_positions =
            topology::delaunayTriangulationVertexPositions( *leaves.back().level_set_tri, level_set_positions );

        process_param( *leaves.back().level_set_tri, tri_positions, thetas, leaves.back() );

        if( log ) std::cout << "FINISHED LEVEL " << ( level_ii + 1 ) << std::endl << std::endl;
    }

    { // target level set
        leaves.push_back( {} );
        const auto target_positions = vertex_positions( bdry );
        leaves.back().reversed_cmap = std::make_unique<topology::ReversedCombinatorialMap>( target );
        const auto& rev_map = *leaves.back().reversed_cmap;
        const auto rev_positions = reversedVertexPositions( rev_map, target_positions );
        const auto face_ids_of_edge = [&]( const topology::Edge& e ) {
            return map_face_ids( bdry.toUnderlyingCell(
                topology::Face( phi( bdry, 2, rev_map.toUnderlyingCell( e ).dart() ).value() ) ) );
        };
        const std::map<topology::Vertex, double> thetas =
            reparam::thetaValues( rev_map, rev_positions, face_ids_of_edge, intersections.back() );

        process_param( rev_map, rev_positions, thetas, leaves.back() );
    }
    if( log ) std::cout << "FINISHED TARGET\n\n";

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
}

void testLevelSetBasedTracing( const SweepInput& sweep_input, const size_t n_levels, const bool log, const std::optional<std::string> output_prefix = std::nullopt )
{
    const std::vector<double> level_set_values = [&]( const size_t n_levels ) {
        std::vector<double> out;
        out.reserve( n_levels );
        const double diff = 1.0 / ( n_levels - 1 );
        for( size_t i = 0; i < n_levels; i++ )
        {
            out.push_back( i * diff );
        }
        return out;
    }( n_levels );

    testLevelSetBasedTracing( sweep_input, level_set_values, log, output_prefix );
}

TEST_CASE( "Level set parameterization of left ventricle" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/left_ventricle.inp", "Surface3", "Surface2" );

    const auto linspace = []( const double left_val, const double right_val, const size_t n_levels ) {
        std::vector<double> out;
        out.reserve( n_levels );
        const double diff = ( right_val - left_val ) / ( n_levels - 1 );
        for( size_t i = 0; i < n_levels; i++ )
        {
            out.push_back( left_val + i * diff );
        }
        return out;
    };

    const auto concatenate = []( const std::vector<double>& first, const std::vector<double>& second ) {
        std::vector<double> out;
        out.reserve( first.size() + second.size() );

        out.insert( out.end(), first.begin(), first.end() );
        out.insert( out.end(), second.begin(), second.end() );

        return out;
    };

    // level set every 0.02 until 0.75, then ~30 levels between 0.75 and 0.85, then every 0.2 until 1.0.
    const std::vector<double> level_set_values =
        concatenate(
            concatenate( linspace( 0, 0.78, 40 ), linspace( 0.8, 0.819, 20 ) ),
            concatenate( concatenate( linspace( 0.81925, 0.81975, 3 ), linspace( 0.82, 0.85, 31 ) ), linspace( 0.86, 1.0, 13 ) ) );

    const bool log_progress = false;
    const std::optional<std::string> output_prefix = std::nullopt;//{ "ventricle" };
    testLevelSetBasedTracing( sweep_input, level_set_values, log_progress, output_prefix );
}

TEST_CASE( "Level set parameterization of hook" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/hook.inp", "Surface12", "Surface10" );

    const bool log_progress = false;
    const std::optional<std::string> output_prefix = std::nullopt;//{ "hook" };
    testLevelSetBasedTracing( sweep_input, 30, log_progress, output_prefix );
}

TEST_CASE( "Level set parameterization of cube" )
{
    const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();

    const bool log_progress = false;
    testLevelSetBasedTracing( sweep_input, 3, log_progress );
}