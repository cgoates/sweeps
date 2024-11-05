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
#include <random>

#include <fstream>
#include <sstream>

std::function<bool( const topology::Vertex& )> testEqualVertices( const topology::IndexingFunc& vert_ids,
                                                                  const topology::Vertex& end_v )
{
    return [&]( const topology::Vertex& test_v ) { return vert_ids( test_v ) == vert_ids( end_v ); };
}

std::function<bool( const topology::Vertex& )> testFoundOtherBoundary( const topology::CombinatorialMap& cmap,
                                                                       const topology::Vertex& v_on_bdry )
{
    auto d = v_on_bdry.dart();
    if( not onBoundary( cmap, d ) )
    {
        bool temp = iterateDartsOfCell( cmap, v_on_bdry, [&]( const auto& a ) {
            d = a;
            return not onBoundary( cmap, a );
        } );
        std::cout << "Did I find the other edge? " << not temp << std::endl;
    }

    const topology::CombinatorialMapBoundary bdry( cmap, { d } );
    const auto bdry_ids = indexingOrError( bdry, 0 );

    std::set<size_t> stop_ids;

    auto curr_d = d;
    do
    {
        stop_ids.insert( bdry_ids( topology::Vertex( curr_d ) ) );
        curr_d = phi( bdry, 1, curr_d ).value();
    } while( curr_d != d );

    const auto v_ids = indexingOrError( cmap, 0 );
    return [stop_ids, v_ids]( const topology::Vertex& v ) { return stop_ids.count( v_ids( v ) ) > 0; };
}

constexpr bool interior_only = true;

reparam::FoliationLeaf
    leafFromLevelSetWithCuts( const std::shared_ptr<const topology::CombinatorialMap>& cmap_ptr,
                              const auto& positions,
                              const std::vector<std::pair<topology::Vertex, topology::Vertex>>& cut_vertices )
{
    if( cut_vertices.size() == 0 ) throw std::invalid_argument( "Must be at least one cut" );
    const auto& cmap = *cmap_ptr;
    const auto cmap_vert_ids = indexingOrError( cmap, 0 );
    std::set<topology::Cell> cuts;
    std::vector<std::pair<topology::Vertex, topology::Vertex>> updated_cut_vertices;
    size_t i = 0;
    for( const auto& pr : cut_vertices )
    {
        const auto cut = topology::shortestPath( cmap, positions, pr.first, testFoundOtherBoundary( cmap, pr.second ), interior_only );
        std::cout << "Cut size: " << cut.size() << std::endl;
        io::outputEdges( cmap, positions, cut, "level_set_cut_" + std::to_string( i++ ) + ".vtu" );
        cuts.insert( cut.begin(), cut.end() );
        updated_cut_vertices.push_back( { pr.first, topology::Vertex( phi( cmap, 1, cut.back().dart() ).value() ) } );
    }

    const auto cut_cmap = std::make_shared<const topology::CutCombinatorialMap>( cmap, cuts );
    const auto cut_vert_ids = indexingOrError( *cut_cmap, 0 );

    const auto is_cut_extremity = [&]( const topology::Vertex& v ) -> bool {
        return std::any_of( updated_cut_vertices.begin(), updated_cut_vertices.end(), [&]( const auto& pr ) {
            return testEqualVertices( cmap_vert_ids, pr.first )( v ) or
                   testEqualVertices( cmap_vert_ids, pr.second )( v );
        } );
    };

    const size_t n_cuts = cut_vertices.size();

    const auto is_start_vert = [&]( const topology::Vertex& v ) {
        const size_t vert_id = cmap_vert_ids( cut_vertices.front().first );
        return cmap_vert_ids( v ) == vert_id and onBoundary( cmap, phi( cmap, -1, v.dart() ).value() );
    };

    const std::map<topology::Vertex, Eigen::Vector2d> bdry_constraints =
        reparam::boundaryConstraints( *cut_cmap, positions, n_cuts, is_cut_extremity, is_start_vert );

    // FIXME: This should be what we return from boundaryConstraints, maybe
    const std::map<size_t, Eigen::Vector2d> constraints_by_id = [&]() {
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
            return constraints_by_id.at( cut_vert_ids( v ) );
        }
        else
            return std::nullopt;
    };

    const Eigen::MatrixX2d tutte = reparam::tutteEmbedding( *cut_cmap, positions, constraints_func, true );

    reparam::FoliationLeaf leaf{};
    leaf.tutte =
        std::make_shared<Eigen::MatrixX2d>( reparam::tutteEmbedding( *cut_cmap, positions, constraints_func, true ) );

    const auto tutte_positions = [cut_vert_ids, tutte = *leaf.tutte]( const topology::Vertex& v ) -> Eigen::Vector2d {
        return ( Eigen::Vector2d() << tutte( cut_vert_ids( v ), 0 ), tutte( cut_vert_ids( v ), 1 ) ).finished();
    };
    const auto cut_atlas = std::make_shared<const param::TriangleParametricAtlas>( cut_cmap );
    const auto atlas = std::make_shared<const param::TriangleParametricAtlas>( cmap_ptr );
    leaf.tutte_mapping = std::make_shared<const mapping::TriangleMeshMapping>( cut_atlas, tutte_positions, 2 );
    leaf.space_mapping = std::make_shared<const mapping::TriangleMeshMapping>( atlas, positions, 3 );
    return leaf;
}

// outer vector is by level set, inner vector is by tunnel loop
std::vector<std::vector<std::pair<topology::Edge, topology::Edge>>>
    tunnelLoopIntersections( const topology::CombinatorialMapBoundary& bdry,
                             const topology::CombinatorialMapRestriction& sides,
                             const std::vector<std::array<topology::Vertex, 4>>& tunnel_loop_vs,
                             const auto& level_set_values,
                             const auto& bdry_positions,
                             const auto& harmonic_func,
                             const std::optional<std::string>& filename = std::nullopt )
{
    const auto side_vertex_ids = indexingOrError( sides, 0 );
    std::vector<std::vector<std::pair<topology::Edge, topology::Edge>>> out(
        level_set_values.size(), std::vector<std::pair<topology::Edge, topology::Edge>>( tunnel_loop_vs.size() ) );

    for( size_t loop_ii = 0; loop_ii < tunnel_loop_vs.size(); loop_ii++ )
    {
        const std::array<topology::Vertex, 4>& vs = tunnel_loop_vs.at( loop_ii );
        const auto path0 = topology::shortestPath(
            sides, bdry_positions, vs.at( 0 ), testEqualVertices( side_vertex_ids, vs.at( 1 ) ) );
        const auto path2 = topology::shortestPath(
            sides, bdry_positions, vs.at( 2 ), testEqualVertices( side_vertex_ids, vs.at( 3 ) ) );

        if( true or filename )
        {
            const auto bdry_vertex_ids = indexingOrError( bdry, 0 );
            const auto path1 = topology::shortestPath(
                bdry, bdry_positions, vs.at( 1 ), testEqualVertices( bdry_vertex_ids, vs.at( 2 ) ) );
            const auto path3 = topology::shortestPath(
                bdry, bdry_positions, vs.at( 3 ), testEqualVertices( bdry_vertex_ids, vs.at( 0 ) ) );
            const auto tunnel_loop =
                util::concatenate( util::concatenate( path0, path1 ), util::concatenate( path2, path3 ) );
            io::outputEdges( bdry, bdry_positions, tunnel_loop, "handle_loop_" + std::to_string( loop_ii ) + ".vtu" );
        }

        std::vector<SmallVector<topology::Edge, 2>> intersections( level_set_values.size(),
                                                                   SmallVector<topology::Edge, 2>() );

        const auto add_any_intersections = [&]( const topology::Edge& e ) {
            const std::pair<double, double> bounds = {
                harmonic_func( bdry.toUnderlyingCell( topology::Vertex( e.dart() ) ) ),
                harmonic_func( bdry.toUnderlyingCell( topology::Vertex( phi( bdry, 1, e.dart() ).value() ) ) ) };
            for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
            {
                if( ( bounds.first < level_set_values.at( level_ii ) and
                      bounds.second >= level_set_values.at( level_ii ) ) or
                    ( bounds.first >= level_set_values.at( level_ii ) and
                      bounds.second < level_set_values.at( level_ii ) ) )
                {
                    intersections.at( level_ii ).push_back( e );
                }
            }
        };

        for( const topology::Edge& e : path0 ) add_any_intersections( e );
        for( const topology::Edge& e : path2 ) add_any_intersections( e );

        for( size_t level_ii = 1; level_ii < intersections.size() - 1; level_ii++ )
            out.at( level_ii - 1 ).at( loop_ii ) = { intersections.at( level_ii ).at( 0 ),
                                                     intersections.at( level_ii ).at( 1 ) };
    }

    return out;
}

void nonDiskFoliations( const SweepInput& sweep,
                        const std::vector<double>& level_set_values,
                        const std::vector<std::array<size_t, 4>>& cut_v_ids,
                        const std::function<void( const std::vector<reparam::FoliationLeaf>& )>& callback )
{
    const topology::TetMeshCombinatorialMap map( sweep.mesh );
    const topology::CombinatorialMapBoundary bdry( map );

    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = reparam::sweepEmbedding( map, sweep.zero_bcs, sweep.one_bcs, normals );

    const auto bdry_vertex_ids = indexingOrError( bdry, 0 );
    const auto map_face_ids = indexingOrError( map, 2 );

    const auto keep_face_base = [&]( const topology::Face& f ) {
        return sweep.zero_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
               sweep.zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
               sweep.zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };

    const auto keep_face_target = [&]( const topology::Face& f ) {
        return sweep.one_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
               sweep.one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
               sweep.one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
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

    const auto side_vertex_ids = indexingOrError( sides, 0 );

    const std::vector<std::array<topology::Vertex, 4>> vs = [&]() {
        std::vector<std::array<topology::Vertex, 4>> out( cut_v_ids.size() );
        iterateCellsWhile( sides, 0, [&]( const topology::Vertex& v ) {
            const size_t id = side_vertex_ids( v );
            for( size_t loop_ii = 0; loop_ii < cut_v_ids.size(); loop_ii++ )
            {
                const auto it = std::find( cut_v_ids.at( loop_ii ).begin(), cut_v_ids.at( loop_ii ).end(), id );
                if( it != cut_v_ids.at( loop_ii ).end() )
                {
                    out.at( loop_ii ).at( std::distance( cut_v_ids.at( loop_ii ).begin(), it ) ) = v;
                }
            }
            return true;
        } );

        return out;
    }();

    const auto harmonic_func = [vol_vertex_ids = indexingOrError( map, 0 ), &sol]( const topology::Vertex& v ) {
        return sol( vol_vertex_ids( v ) );
    };

    const std::vector<std::vector<std::pair<topology::Edge, topology::Edge>>> tunnel_loop_intersections =
        tunnelLoopIntersections( bdry, sides, vs, level_set_values, vertex_positions( bdry ), harmonic_func );

    using namespace topology;
    std::vector<reparam::FoliationLeaf> leaves;

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
        std::vector<std::pair<topology::Vertex, topology::Vertex>> cut_ends;
        cut_ends.reserve( vs.size() );
        for( const auto& loop_vs : vs )
        {
            cut_ends.push_back( { find_on_base( loop_vs.at( 0 ) ), find_on_base( loop_vs.at( 3 ) ) } );
        }
        const auto base_positions = vertex_positions( bdry );
        leaves.push_back( leafFromLevelSetWithCuts( base, base_positions, cut_ends ) );
    }

    for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
    { // Midway level set
        const auto level_set = std::make_shared<const topology::LevelSetCMap>( map, harmonic_func, level_set_values[level_ii] );
        const auto v_pos = vertex_positions( map );
        const auto level_set_positions = topology::levelSetVertexPositions( *level_set, v_pos );
        const auto level_set_tri = std::make_shared<const topology::DelaunayTriangulation>( level_set, level_set_positions );
        const auto tri_positions =
            topology::delaunayTriangulationVertexPositions( *level_set_tri, level_set_positions );

        std::vector<std::pair<topology::Vertex, topology::Vertex>> cut_ends;
        cut_ends.reserve( vs.size() );
        for( const auto& loop_ends : tunnel_loop_intersections.at( level_ii - 1 ) )
        {
            cut_ends.push_back(
                { topology::Vertex( loop_ends.first.dart() ), topology::Vertex( loop_ends.second.dart() ) } );
        }

        leaves.push_back( leafFromLevelSetWithCuts( level_set_tri, tri_positions, cut_ends ) );
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
        std::vector<std::pair<topology::Vertex, topology::Vertex>> cut_ends;
        cut_ends.reserve( vs.size() );
        for( const auto& loop_vs : vs )
        {
            cut_ends.push_back( { rev_map->fromUnderlyingCell( find_on_target( loop_vs.at( 1 ) ) ),
                                  rev_map->fromUnderlyingCell( find_on_target( loop_vs.at( 2 ) ) ) } );
        }
        leaves.push_back( leafFromLevelSetWithCuts( rev_map, rev_positions, cut_ends ) );
    }

    callback( leaves );
}

void macaroniFoliation( const std::vector<double>& level_set_values,
                        const std::function<void( const std::vector<reparam::FoliationLeaf>& )>& callback )
{
    const SweepInput sweep = io::loadINPFile( SRC_HOME "/test/data/macaroni_coarse.inp", "Surface3", "Surface4" );

    const std::array<size_t, 4> cut_v_ids{ 97, 185, 69, 13 }; // base, target, target, base

    nonDiskFoliations( sweep, level_set_values, { cut_v_ids }, callback );
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