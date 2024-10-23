#include <catch2/catch_test_macros.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <TriangleParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <TriangleMeshMapping.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>
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
#include <GeometricMapping.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <LeastSquaresFitting.hpp>
#include <MultiPatchSplineSpace.hpp>

void macaroniFoliation( const std::vector<double>& level_set_values,
                        const std::function<void( const std::vector<reparam::FoliationLeaf>& )>& callback )
{
    const SweepInput sweep = io::loadINPFile( SRC_HOME "/test/data/macaroni_coarse.inp", "Surface3", "Surface4" );
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

    std::vector<reparam::FoliationLeaf> leaves;

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



        leaves.push_back({});
        leaves.back().tutte = std::make_shared<Eigen::MatrixX2d>( reparam::tutteEmbedding( *cut_cmap, positions, constraints_func, true ) );

        const auto tutte_positions = [cut_vert_ids,tutte=*leaves.back().tutte]( const topology::Vertex& v ) -> Eigen::Vector2d {
            return ( Eigen::Vector2d() << tutte( cut_vert_ids( v ), 0 ), tutte( cut_vert_ids( v ), 1 ) ).finished();
        };
        const auto cut_atlas = std::make_shared<const param::TriangleParametricAtlas>( cut_cmap );
        const auto atlas = std::make_shared<const param::TriangleParametricAtlas>( cmap_ptr );
        leaves.back().tutte_mapping = std::make_shared<const mapping::TriangleMeshMapping>( cut_atlas, tutte_positions, 2 );
        leaves.back().space_mapping = std::make_shared<const mapping::TriangleMeshMapping>( atlas, positions, 3 );
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

    callback( leaves );
}

TEST_CASE( "Level set parameterization of the macaroni" )
{
    const size_t n_levels = 25;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    SimplicialComplex level_sets;
    const std::vector<Eigen::Vector2d> square_points{
        { 0.00000000e+00, -9.00000000e-01 },  { -2.25000000e-01, -6.75000000e-01 },
        { -4.50000000e-01, -4.50000000e-01 }, { -6.75000000e-01, -2.25000000e-01 },
        { -9.00000000e-01, 0.00000000e+00 },  { 2.25000000e-01, -6.75000000e-01 },
        { 0.00000000e+00, -4.50000000e-01 },  { -2.25000000e-01, -2.25000000e-01 },
        { -4.50000000e-01, 0.00000000e+00 },  { -6.75000000e-01, 2.25000000e-01 },
        { 4.50000000e-01, -4.50000000e-01 },  { 2.25000000e-01, -2.25000000e-01 },
        { 0.00000000e+00, 0.00000000e+00 },   { -2.25000000e-01, 2.25000000e-01 },
        { -4.50000000e-01, 4.50000000e-01 },  { 6.75000000e-01, -2.25000000e-01 },
        { 4.50000000e-01, 2.77555756e-17 },   { 2.25000000e-01, 2.25000000e-01 },
        { 0.00000000e+00, 4.50000000e-01 },   { -2.25000000e-01, 6.75000000e-01 },
        { 9.00000000e-01, 0.00000000e+00 },   { 6.75000000e-01, 2.25000000e-01 },
        { 4.50000000e-01, 4.50000000e-01 },   { 2.25000000e-01, 6.75000000e-01 },
        { 0.00000000e+00, 9.00000000e-01 } };

    std::vector<std::vector<Eigen::Vector3d>> mapped_points( square_points.size() );

    macaroniFoliation( level_set_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
        for( size_t level_ii = 0; level_ii < leaves.size(); level_ii++ )
        {
            const auto& leaf = leaves.at( level_ii );
            const auto& square_mapping = leaf.tutte_mapping;
            const auto& cmap = square_mapping->parametricAtlas().cmap();
            const auto& positions = leaf.space_mapping->vertPositions();
            for( size_t point_ii = 0; point_ii < square_points.size(); point_ii++ )
            {
                const auto param_pt = square_mapping->maybeInverse( square_points.at( point_ii ) );
                CHECK( param_pt.has_value() );
                if( not param_pt.has_value() )
                {
                    std::cerr << "NO VALUE level: " << level_ii << " pt: " << square_points.at( point_ii ) << std::endl;
                    break;
                }
                const auto space_pt = leaf.space_mapping->evaluate( param_pt.value().first, param_pt.value().second );
                mapped_points.at( point_ii ).push_back( space_pt );
            }

            iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
                addTriangleNoDuplicateChecking( level_sets, triangleOfFace<3>( cmap, positions, f ) );
                return true;
            } );
        }
    } );

    io::VTKOutputObject output( level_sets );
    io::outputSimplicialFieldToVTK( output, "macaroni_level_sets.vtu" );

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

Eigen::MatrixXd fitToLeaves(
    eval::SplineSpaceEvaluator& evaler,
    const std::vector<double>& level_set_values,
    const size_t n_points,
    const std::function<void(
        const size_t,
        const std::function<Eigen::Vector3d( const topology::Cell&, const param::ParentPoint&, const Eigen::Vector2d& )>& )>&
        leaf_point_iterator )
{
    const param::ParentDomain pd_3d = param::cubeDomain( 3 );

    // If not covering the whole domain, add source and target surface
    std::vector<double> level_set_values2;
    const bool add_front = not util::equals( level_set_values.front(), 0.0, 1e-9 );
    const bool add_back = not util::equals( level_set_values.back(), 1.0, 1e-9 );
    if( add_front ) level_set_values2.push_back( 0 );
    level_set_values2.insert( level_set_values2.end(), level_set_values.begin(), level_set_values.end() );
    if( add_back ) level_set_values2.push_back( 1.0 );

    Eigen::MatrixXd fit_cpts;

    macaroniFoliation( level_set_values2, [&]( const std::vector<reparam::FoliationLeaf>& leaves2 ) {
        const std::vector<reparam::FoliationLeaf> leaves( add_front ? std::next( leaves2.begin() ) : leaves2.begin(), add_back ? std::prev( leaves2.end() ) : leaves2.end() );
        std::cout << "About to fit\n";

        fit_cpts = fitting::leastSquaresFitting(
            evaler,
            n_points,
            3,
            [&]( const std::function<void( const topology::Cell&, const param::ParentPoint&, const Eigen::VectorXd& )>&
                    add_least_squares_point ) {
                for( size_t leaf_ii = 0; leaf_ii < leaves.size(); leaf_ii++ )
                {
                    leaf_point_iterator( leaf_ii, [&]( const topology::Cell& cell, const param::ParentPoint& vol_ppt, const Eigen::Vector2d& circle_pt ) {
                        const auto param_pt = leaves.at( leaf_ii ).tutte_mapping->maybeInverse( circle_pt );
                        CHECK( param_pt.has_value() );
                        if( not param_pt.has_value() )
                        {
                            std::cerr << "NO VALUE" << std::endl;
                            pauseDebugger();
                        }
                        const Eigen::Vector3d field_pt = leaves.at( leaf_ii ).space_mapping->evaluate( param_pt.value().first, param_pt.value().second );

                        add_least_squares_point( cell, vol_ppt, field_pt );
                        return field_pt;
                    } );
                }
            } );
    } );

    return fit_cpts;
}

Eigen::Vector2d toUnitSquare( const topology::TPCombinatorialMap& cmap, const topology::Face& f, const Eigen::Vector2d& pt )
{
    // FIXME: Assumes unit parametric lengths
    const auto [d1, d2, _] = cmap.unflatten( f.dart() );
    const Eigen::Vector2d zero_to_one( ( (double)d1.id() + pt( 0 ) ) / (double)cellCount( cmap.sourceCMap(), 1 ),
                                       ( (double)d2.id() + pt( 1 ) ) / (double)cellCount( cmap.lineCMap(), 1 ) );

    return Eigen::Rotation2D( -1*std::numbers::pi/4 ).toRotationMatrix() * ( sqrt( 2 ) * ( zero_to_one - Eigen::Vector2d::Constant( 0.5 ) ) );
}

void fitToPringlesPeriodicPatch( const std::vector<double>& level_set_values,
                                 const std::string& output_prefix,
                                 const size_t degree,
                                 const basis::KnotVector& kv_u,
                                 const size_t n_elems_s,
                                 const size_t n_elems_t )
{
    using namespace topology;

    const basis::KnotVector kv_s = basis::integerKnotsWithNElems( n_elems_s, degree );
    const basis::KnotVector kv_t = basis::integerKnotsWithNElems( n_elems_t, degree );

    const std::shared_ptr<const basis::TPSplineSpace> TP_ss = std::make_shared<const basis::TPSplineSpace>(
        basis::buildBSpline( { kv_s, kv_t, kv_u }, { degree, degree, degree } ) );
    const auto& vol_cmap = TP_ss->basisComplex().parametricAtlas().cmap();
    const basis::TPSplineSpace& source_ss = static_cast<const basis::TPSplineSpace&>( TP_ss->source() );
    const std::vector<std::pair<topology::Cell, param::ParentPoint>> ppt_u =
        parentPointsOfParamPoints( level_set_values, TP_ss->line().basisComplex().parametricAtlas(), 1e-9 );

    const basis::MultiPatchSplineSpace mp_ss =
        basis::buildMultiPatchSplineSpace( std::vector<std::shared_ptr<const basis::TPSplineSpace>>( 1, TP_ss ),
                                           { { { 0, Dart( 3 * 6 + 1 ) }, { 0, Dart( ( 4 * ( n_elems_s - 1 ) + 1 ) * 6 + 1 ) } } } );

    const SmallVector<double, 2> source_points{ 0.1, 0.9 };
    const param::ParentDomain pd_3d = param::cubeDomain( 3 );

    SimplicialComplex fitting_points;
    eval::SplineSpaceEvaluator evaler( mp_ss, 0 );

    const size_t n_bdry_pts = 8 * level_set_values.size();

    const size_t n_points = cellCount( source_ss.basisComplex().parametricAtlas().cmap(), 2 ) * 4 * 5 * level_set_values.size() + n_bdry_pts;

    Eigen::MatrixXd fit_cpts = fitToLeaves( evaler, level_set_values, n_points, [&]( const size_t leaf_ii, const auto& point_callback ) {
        const size_t patch_ii = 0;
        iterateCellsWhile(
            source_ss.basisComplex().parametricAtlas().cmap(), 2, [&]( const topology::Cell& f ) {
                util::iterateTensorProduct(
                    { source_points.size(), source_points.size() }, [&]( const util::IndexVec& indices ) {
                        const Eigen::Vector3d pt( source_points.at( indices.at( 0 ) ),
                                                source_points.at( indices.at( 1 ) ),
                                                ppt_u.at( leaf_ii ).second.mPoint( 0 ) );
                        const param::ParentPoint vol_ppt(
                            pd_3d,
                            pt,
                            { false,
                            false,
                            false,
                            false,
                            ppt_u.at( leaf_ii ).second.mBaryCoordIsZero.at( 0 ),
                            ppt_u.at( leaf_ii ).second.mBaryCoordIsZero.at( 1 ) } );

                        const topology::Volume patch_cell(
                            vol_cmap.flatten( f.dart(),
                                            ppt_u.at( leaf_ii ).first.dart(),
                                            topology::TPCombinatorialMap::TPDartPos::DartPos0 ) );

                        const topology::Volume cell( mp_ss.basisComplex().parametricAtlas().cmap().toGlobalDart( patch_ii, patch_cell.dart() ) );

                        const Eigen::Vector2d tutte_pt = toUnitSquare( source_ss.basisComplex().parametricAtlas().cmap(), f, pt.head( 2 ) );
                        
                        const Eigen::Vector3d field_pt = point_callback( cell, vol_ppt, tutte_pt );
                        fitting_points.simplices.emplace_back( fitting_points.points.size() );
                        fitting_points.points.push_back( field_pt );
                    } );
                return true;
            } );
    } );

    io::VTKOutputObject fitting_points_output( fitting_points );
    io::outputSimplicialFieldToVTK( fitting_points_output, output_prefix + "fitting_points.vtu" );

    io::outputBezierMeshToVTK( mp_ss, fit_cpts, "fit_to_" + output_prefix + "_multi_patch.vtu" );
}

TEST_CASE( "Spline fit to the macaroni" )
{
    const size_t n_levels = 25;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    const basis::KnotVector kv_u = basis::integerKnotsWithNElems( 10, 2 );
    fitToPringlesPeriodicPatch( level_set_values, "macaroni", 2, kv_u, 10, 2 );
}