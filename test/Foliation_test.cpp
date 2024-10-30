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
#include <Dijkstra.hpp>
#include <CutCombinatorialMap.hpp>

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
                    const auto& param_pt = leaf.tutte_mapping->maybeInverse( circle_pt );
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

TEST_CASE( "Level set parameterization of stanford bunny" )
{
    const SweepInput sweep_input = SweepInputTestCases::bunny();

    using util::linspace;
    using util::concatenate;

    const std::vector<double> level_set_values = concatenate( { linspace( 0, 0.177, 35 ),
                                                                linspace( 0.177, 0.17904, 40 ),
                                                                linspace( 0.17904, 0.17908, 40 ),
                                                                linspace( 0.17908, 0.1796, 10 ),
                                                                linspace( 0.1796, 1.0, 20 ) } );

    const std::optional<std::string> output_prefix = std::nullopt;//{ "bunny" };
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

TEST_CASE( "Level set parameterization of bullet" )
{
    const SweepInput sweep_input = SweepInputTestCases::bullet();

    const std::optional<std::string> output_prefix = std::nullopt;//{ "bullet" };
    testLevelSetBasedTracing( sweep_input, 5, "bullet" );
}

TEST_CASE( "Level set parameterization of flange" )
{
    const size_t n_levels = 25;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    const SweepInput sweep = SweepInputTestCases::flange();
    const auto& zero_bcs = sweep.zero_bcs;
    const auto& one_bcs = sweep.one_bcs;

    const topology::TetMeshCombinatorialMap map( sweep.mesh );
    const topology::CombinatorialMapBoundary bdry( map );

    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = reparam::sweepEmbedding( map, zero_bcs, one_bcs, normals );

    const auto bdry_vertex_ids = indexingOrError( bdry, 0 );
    const auto map_face_ids = indexingOrError( map, 2 );

    const auto keep_face_base = [&]( const topology::Face& f ) {
        return zero_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
            zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
            zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };

    const auto keep_face_target = [&]( const topology::Face& f ) {
        return one_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
            one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
            one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };

    const topology::CombinatorialMapRestriction base( bdry, keep_face_base, true );
    const topology::CombinatorialMapRestriction target( bdry, keep_face_target, true );

    const auto vertex_positions = [&sweep]( const topology::CombinatorialMap& map ) {
        const auto vertex_ids = indexingOrError( map, 0 );
        return [&sweep, vertex_ids]( const topology::Vertex& v ) -> Eigen::Vector3d {
            return sweep.mesh.points.at( vertex_ids( v ) );
        };
    };

    const auto vertex_ids = indexingOrError( map, 0 );

    const auto harmonic_func = [&]( const topology::Vertex& v ) { return sol( vertex_ids( v ) ); };

    SimplicialComplex level_sets;
    const auto process_param = [&]( const topology::CombinatorialMap& cmap, const auto& positions ) {
        iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
            addTriangleNoDuplicateChecking( level_sets, triangleOfFace<3>( cmap, positions, f ) );
            return true;
        } );
    };

    { // Base level set
        const auto base_positions = vertex_positions( bdry );
        process_param( base, base_positions );
    }

    for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
    { // Midway level set
        const auto level_set = std::make_shared<const topology::LevelSetCMap>( map, harmonic_func, level_set_values[level_ii] );
        const auto v_pos = vertex_positions( map );
        const auto level_set_positions = topology::levelSetVertexPositions( *level_set, v_pos );
        const auto level_set_tri =
            std::make_unique<topology::DelaunayTriangulation>( level_set, level_set_positions );
        const auto tri_positions =
            topology::delaunayTriangulationVertexPositions( *level_set_tri, level_set_positions );
        process_param( *level_set_tri, tri_positions );
    }

    { // target level set
        const auto target_positions = vertex_positions( bdry );
        const topology::ReversedCombinatorialMap rev_map( target );
        const auto rev_positions = reversedVertexPositions( rev_map, target_positions );
        process_param( rev_map, rev_positions );
    }

    io::VTKOutputObject output( level_sets );
    io::outputSimplicialFieldToVTK( output, "flange_level_sets.vtu" );
}

TEST_CASE( "Level set parameterization of the macaroni" )
{
    const size_t n_levels = 25;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    const SweepInput sweep = io::loadINPFile( SRC_HOME "/test/data/macaroni_coarse.inp", "Surface3", "Surface4" );// 97, 185, 69, 13
    const auto& zero_bcs = sweep.zero_bcs;
    const auto& one_bcs = sweep.one_bcs;

    const topology::TetMeshCombinatorialMap map( sweep.mesh );
    const topology::CombinatorialMapBoundary bdry( map );

    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = reparam::sweepEmbedding( map, zero_bcs, one_bcs, normals );

    const auto bdry_vertex_ids = indexingOrError( bdry, 0 );
    const auto map_face_ids = indexingOrError( map, 2 );

    const auto keep_face_base = [&]( const topology::Face& f ) {
        return zero_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
            zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
            zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };

    const auto keep_face_target = [&]( const topology::Face& f ) {
        return one_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
            one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
            one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };

    const auto keep_face_sides = [&]( const topology::Face& f ) {
        return not keep_face_base( f ) and not keep_face_target( f );
    };

    const auto base = std::make_shared<const topology::CombinatorialMapRestriction>( bdry, keep_face_base, true );
    const topology::CombinatorialMapRestriction target( bdry, keep_face_target, true );
    const topology::CombinatorialMapRestriction sides( bdry, keep_face_sides );

    const auto vertex_positions = [&sweep]( const topology::CombinatorialMap& map ) {
        const auto vertex_ids = indexingOrError( map, 0 );
        return [&sweep, vertex_ids]( const topology::Vertex& v ) -> Eigen::Vector3d {
            return sweep.mesh.points.at( vertex_ids( v ) );
        };
    };

    const auto test_equal_vertices = []( const topology::IndexingFunc& vert_ids, const topology::Vertex& end_v ) {
        return [&]( const topology::Vertex& test_v ){
            return vert_ids( test_v ) == vert_ids( end_v );
        };
    };

    const auto side_vertex_ids = indexingOrError( sides, 0 );

    const std::array<topology::Vertex, 4> vs = [&]() {
        topology::Vertex v0;
        topology::Vertex v1;
        topology::Vertex v2;
        topology::Vertex v3;
        iterateCellsWhile( sides, 0, [&]( const topology::Vertex& v ) {
            const size_t id = side_vertex_ids( v );
            switch( id )
            {
                case 97:
                    v0 = v;
                    return true;
                case 185:
                    v1 = v;
                    return true;
                case 69:
                    v2 = v;
                    return true;
                case 13:
                    v3 = v;
                    return true;
                default:
                    return true;
            }
        } );
        
        return std::array<topology::Vertex, 4>{ v0, v1, v2, v3 };
    }();
    const auto path0 = topology::shortestPath( bdry, vertex_positions( bdry ), vs.at( 0 ), test_equal_vertices( bdry_vertex_ids, vs.at( 1 ) ) );
    const auto path1 = topology::shortestPath( bdry, vertex_positions( bdry ), vs.at( 1 ), test_equal_vertices( bdry_vertex_ids, vs.at( 2 ) ) );
    const auto path2 = topology::shortestPath( sides, vertex_positions( sides ), vs.at( 2 ), test_equal_vertices( side_vertex_ids, vs.at( 3 ) ) );
    const auto path3 = topology::shortestPath( bdry, vertex_positions( bdry ), vs.at( 3 ), test_equal_vertices( bdry_vertex_ids, vs.at( 0 ) ) );
    const auto tunnel_loop = util::concatenate( util::concatenate( path0, path1 ), util::concatenate( path2, path3 ) );

    [&]( const std::vector<topology::Edge>& edges, const std::string& filename ) {
        SimplicialComplex path;

        const auto bdry_positions = vertex_positions( bdry );

        path.points.push_back( bdry_positions( topology::Vertex( edges.front().dart() ) ) );
        for( const auto& edge : edges )
        {
            const size_t temp = path.points.size();
            path.points.push_back( bdry_positions( topology::Vertex( phi( bdry, 1, edge.dart() ).value() ) ) );
            path.simplices.emplace_back( temp - 1, temp );
        }

        io::VTKOutputObject output( path );
        io::outputSimplicialFieldToVTK( output, filename );
    }( tunnel_loop, "macaroni_tunnel_loop.vtu" );

    const auto vol_vertex_ids = indexingOrError( map, 0 );

    const auto harmonic_func = [&]( const topology::Vertex& v ) { return sol( vol_vertex_ids( v ) ); };

    const std::vector<std::pair<topology::Edge, topology::Edge>> tunnel_loop_intersections = [&](){
        std::vector<SmallVector<topology::Edge, 2>> intersections( level_set_values.size(), SmallVector<topology::Edge, 2>() );
        for( const topology::Edge& e : tunnel_loop )
        {
            const std::pair<double, double> bounds = {
                sol( bdry_vertex_ids( topology::Vertex( e.dart() ) ) ),
                sol( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, e.dart() ).value() ) ) ) };
            for( size_t i = 1; i < level_set_values.size() - 1; i++ )
            {
                if( ( bounds.first < level_set_values.at( i ) and bounds.second >= level_set_values.at( i ) ) or
                    ( bounds.first >= level_set_values.at( i ) and bounds.second < level_set_values.at( i ) ) )
                {
                    intersections.at( i ).push_back( e );
                }
            }
        }
        std::vector<std::pair<topology::Edge, topology::Edge>> out;
        out.reserve( intersections.size() - 2 );
        for( size_t i = 1; i < intersections.size() - 1; i++ )
            out.emplace_back( intersections.at( i ).at( 0 ), intersections.at( i ).at( 1 ) );

        return out;
    }();

    SimplicialComplex cuts;
    SimplicialComplex level_sets;
    const std::vector<Eigen::Vector2d> square_points{
       { 0.00000000e+00, -9.00000000e-01},
       {-2.25000000e-01, -6.75000000e-01},
       {-4.50000000e-01, -4.50000000e-01},
       {-6.75000000e-01, -2.25000000e-01},
       {-9.00000000e-01,  0.00000000e+00},
       { 2.25000000e-01, -6.75000000e-01},
       { 0.00000000e+00, -4.50000000e-01},
       {-2.25000000e-01, -2.25000000e-01},
       {-4.50000000e-01,  0.00000000e+00},
       {-6.75000000e-01,  2.25000000e-01},
       { 4.50000000e-01, -4.50000000e-01},
       { 2.25000000e-01, -2.25000000e-01},
       { 0.00000000e+00,  0.00000000e+00},
       {-2.25000000e-01,  2.25000000e-01},
       {-4.50000000e-01,  4.50000000e-01},
       { 6.75000000e-01, -2.25000000e-01},
       { 4.50000000e-01,  2.77555756e-17},
       { 2.25000000e-01,  2.25000000e-01},
       { 0.00000000e+00,  4.50000000e-01},
       {-2.25000000e-01,  6.75000000e-01},
       { 9.00000000e-01,  0.00000000e+00},
       { 6.75000000e-01,  2.25000000e-01},
       { 4.50000000e-01,  4.50000000e-01},
       { 2.25000000e-01,  6.75000000e-01},
       { 0.00000000e+00,  9.00000000e-01}
    };

    std::vector<std::vector<Eigen::Vector3d>> mapped_points( square_points.size() );

    using namespace topology;
    size_t level_ii = 0;
    const auto process_param = [&]( const std::shared_ptr<const topology::CombinatorialMap>& cmap_ptr,
                                    const auto& positions,
                                    const Vertex& v0,
                                    const Vertex& v1 ) {
        const auto& cmap = *cmap_ptr;
        const auto cmap_vert_ids = indexingOrError( cmap, 0 );
        const auto cut = topology::shortestPath( cmap, positions, v0, test_equal_vertices( cmap_vert_ids, v1 ) );

        const auto cut_cmap = std::make_shared<const topology::CutCombinatorialMap>(
            cmap, std::set<topology::Cell>( cut.begin(), cut.end() ) );
        const auto cut_vert_ids = indexingOrError( *cut_cmap, 0 );

        const auto is_cut_extremity = [&]( const topology::Vertex& v ) -> bool {
            return test_equal_vertices( cmap_vert_ids, v0 )( v ) or test_equal_vertices( cmap_vert_ids, v1 )( v );
        };

        const size_t n_cuts = 1;

        const auto is_start_vert = [&]( const topology::Vertex& v ) {
            const size_t vert_id = cmap_vert_ids( v0 );
            return cmap_vert_ids( v ) == vert_id and onBoundary( cmap, phi( cmap, -1, v.dart() ).value() );
        };

        const std::map<topology::Vertex, Eigen::Vector2d> bdry_constraints =
            reparam::boundaryConstraints( *cut_cmap, positions, n_cuts, is_cut_extremity, is_start_vert );

        const std::map<size_t, Eigen::Vector2d> thetas_by_id = [&]() {
            std::map<size_t, Eigen::Vector2d> out;
            for( const auto& pr : bdry_constraints )
            {
                out.insert( { cut_vert_ids( pr.first ), pr.second } );
            }
            return out;
        }();

        const auto constraints_func = [&]( const topology::Vertex& v ) -> std::optional<Eigen::Vector2d> {
            if( boundaryAdjacent( *cut_cmap, v ) )
            {
                return thetas_by_id.at( cut_vert_ids( v ) );
            }
            else
                return std::nullopt;
        };

        const Eigen::MatrixX2d tutte = reparam::tutteEmbedding( *cut_cmap, positions, constraints_func, true );

        const auto tutte_positions = [&]( const topology::Vertex& v ) -> Eigen::Vector2d {
            return ( Eigen::Vector2d() << tutte( cut_vert_ids( v ), 0 ), tutte( cut_vert_ids( v ), 1 ) ).finished();
        };

        const auto cut_atlas = std::make_shared<const param::TriangleParametricAtlas>( cut_cmap );
        const auto atlas = std::make_shared<const param::TriangleParametricAtlas>( cmap_ptr );
        const auto square_mapping = std::make_shared<const mapping::TriangleMeshMapping>( cut_atlas, tutte_positions, 2 );
        const auto space_mapping = std::make_shared<const mapping::TriangleMeshMapping>( atlas, positions, 3 );

        { //OUTPUT
            for( size_t point_ii = 0; point_ii < square_points.size(); point_ii++ )
            {
                const auto param_pt = square_mapping->maybeInverse( square_points.at( point_ii ) );
                CHECK( param_pt.has_value() );
                if( not param_pt.has_value() )
                {
                    std::cerr << "NO VALUE level: " << level_ii << " pt: " << square_points.at( point_ii ) << std::endl;
                    break;
                }
                const auto space_pt = space_mapping->evaluate( param_pt.value().first, param_pt.value().second );
                mapped_points.at( point_ii ).push_back( space_pt );
            }

            iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
                addTriangleNoDuplicateChecking( level_sets, triangleOfFace<3>( cmap, positions, f ) );
                return true;
            } );

            cuts.points.push_back( positions( topology::Vertex( cut.front().dart() ) ) );
            for( const auto& edge : cut )
            {
                const size_t temp = cuts.points.size();
                cuts.points.push_back( positions( topology::Vertex( phi( cmap, 1, edge.dart() ).value() ) ) );
                cuts.simplices.emplace_back( temp - 1, temp );
            }
        }
        level_ii++;
    };

    { // Base level set
        const auto find_on_base = [&]( const Vertex& v ) {
            std::optional<Vertex> out;
            iterateDartsOfCell( bdry, v, [&]( const Dart& d ) {
                if( keep_face_base( d ) )
                {
                    out.emplace( d );
                    return false;
                }
                return true;
            } );
            return out.value();
        };
        const Vertex v0 = find_on_base( vs.at( 0 ) );
        const Vertex v1 = find_on_base( vs.at( 3 ) );
        const auto base_positions = vertex_positions( bdry );
        process_param( base, base_positions, v0, v1 );
    }

    for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
    { // Midway level set
        const auto level_set = std::make_shared<const topology::LevelSetCMap>( map, harmonic_func, level_set_values[level_ii] );
        const auto v_pos = vertex_positions( map );
        const auto level_set_positions = topology::levelSetVertexPositions( *level_set, v_pos );
        const auto level_set_tri = std::make_shared<const topology::DelaunayTriangulation>( level_set, level_set_positions );
        const auto tri_positions =
            topology::delaunayTriangulationVertexPositions( *level_set_tri, level_set_positions );

        const topology::Vertex v0( tunnel_loop_intersections.at( level_ii - 1 ).first.dart() );
        const topology::Vertex v1( tunnel_loop_intersections.at( level_ii - 1 ).second.dart() );
        process_param( level_set_tri, tri_positions, v0, v1 );
    }

    { // target level set
        const auto find_on_target = [&]( const Vertex& v ) {
            std::optional<Vertex> out;
            iterateDartsOfCell( bdry, v, [&]( const Dart& d ) {
                if( keep_face_target( d ) )
                {
                    out.emplace( d );
                    return false;
                }
                return true;
            } );
            return out.value();
        };
        
        const auto target_positions = vertex_positions( bdry );
        const auto rev_map = std::make_shared<const topology::ReversedCombinatorialMap>( target );
        const auto rev_positions = reversedVertexPositions( *rev_map, target_positions );
        const Vertex v0 = rev_map->fromUnderlyingCell( find_on_target( vs.at( 1 ) ) );
        const Vertex v1 = rev_map->fromUnderlyingCell( find_on_target( vs.at( 2 ) ) );
        process_param( rev_map, rev_positions, v0, v1 );
    }

    io::VTKOutputObject output( level_sets );
    io::outputSimplicialFieldToVTK( output, "macaroni_level_sets.vtu" );

    io::VTKOutputObject output2( cuts );
    io::outputSimplicialFieldToVTK( output2, "macaroni_level_set_cuts.vtu" );

    SimplicialComplex traces;
    for( const auto& line : mapped_points )
    {
        traces.points.push_back( line.front() );
        for( size_t i = 1; i < line.size(); i++ )
        {
            traces.points.push_back( line.at( i ) );
            traces.simplices.emplace_back( traces.points.size() - 1, traces.points.size() - 2 );
        }
    }
    io::VTKOutputObject output4( traces );
    io::outputSimplicialFieldToVTK( output4, "macaroni_traces.vtu" );
}