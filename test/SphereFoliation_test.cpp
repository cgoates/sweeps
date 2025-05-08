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
#include <GlobalCellMarker.hpp>

#include <fstream>
#include <sstream>

#include <Eigen/Dense>

std::function<bool( const topology::Vertex& )> testEqualVertices( const topology::IndexingFunc& vert_ids,
                                                                  const topology::Vertex& end_v )
{
    return [&]( const topology::Vertex& test_v ) { return vert_ids( test_v ) == vert_ids( end_v ); };
}

static int j = 0;
reparam::FoliationLeaf leafFromSphereLevelSet( const std::shared_ptr<const topology::CombinatorialMap>& cmap_ptr,
                                               const VertexPositionsFunc& positions,
                                               const std::array<topology::Vertex, 3>& cut_vertices )
{
    const auto& cmap = *cmap_ptr;
    const auto cmap_vert_ids = indexingOrError( cmap, 0 );
    std::set<topology::Cell> cuts;
    size_t i = 0;

    const auto cut1 = topology::shortestPath( cmap, positions, cut_vertices.at( 0 ), testEqualVertices( cmap_vert_ids, cut_vertices.at( 1 ) ) );
    topology::GlobalCellMarker cut_marker( cmap, 1 );
    for( const auto& e : cut1 ) cut_marker.mark( cmap, e );
    const auto cut2 = topology::shortestPath( cmap, [&]( const topology::Edge& e ){
        return cut_marker.isMarked( e ) ? std::numeric_limits<double>::max() : edgeLength( cmap, positions, e );
    }, cut_vertices.at( 1 ), testEqualVertices( cmap_vert_ids, cut_vertices.at( 2 ) ) );

    io::outputEdges( cmap, positions, util::concatenate( cut1, cut2 ), "level_set_cut_" + std::to_string( j ) + ".vtu" );
    j++;
    cuts.insert( cut1.begin(), cut1.end() );
    cuts.insert( cut2.begin(), cut2.end() );

    const auto cut_cmap = std::make_shared<const topology::CutCombinatorialMap>( cmap, cuts );
    const auto cut_vert_ids = indexingOrError( *cut_cmap, 0 );

    const auto is_cut_extremity = [&]( const topology::Vertex& v ) -> bool {
        return std::any_of( cut_vertices.begin(), cut_vertices.end(), [&]( const auto& pr ) {
            return testEqualVertices( cmap_vert_ids, pr )( v );
        } );
    };

    const size_t n_cuts = 1;

    const auto is_start_vert = [&]( const topology::Vertex& v ) {
        const size_t vert_id = cmap_vert_ids( cut_vertices.front() );
        return cmap_vert_ids( v ) == vert_id;
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

// outer vector is by level set
std::vector<std::array<topology::Edge, 3>>
    traceIntersections( const topology::CombinatorialMap& map,
                        const std::array<std::array<topology::Vertex, 3>, 2>& tunnel_loop_vs,
                        const std::vector<double>& level_set_values,
                        const VertexPositionsFunc& map_positions,
                        const std::function<double( const topology::Vertex& )>& harmonic_func,
                        const std::optional<std::string>& filename = std::nullopt )
{
    const auto vertex_ids = indexingOrError( map, 0 );
    const auto path0 = topology::shortestPath(
        map, map_positions, tunnel_loop_vs.at( 0 ).at( 0 ), testEqualVertices( vertex_ids, tunnel_loop_vs.at( 1 ).at( 0 ) ) );
    const auto path1 = topology::shortestPath(
        map, map_positions, tunnel_loop_vs.at( 0 ).at( 1 ), testEqualVertices( vertex_ids, tunnel_loop_vs.at( 1 ).at( 1 ) ) );
    const auto path2 = topology::shortestPath(
        map, map_positions, tunnel_loop_vs.at( 0 ).at( 2 ), testEqualVertices( vertex_ids, tunnel_loop_vs.at( 1 ).at( 2 ) ) );

    if( true or filename )
    {
        io::outputEdges( map, map_positions, path0, "sphere_cut_lines0.vtu" );
        io::outputEdges( map, map_positions, path1, "sphere_cut_lines1.vtu" );
        io::outputEdges( map, map_positions, path2, "sphere_cut_lines2.vtu" );
    }

    std::vector<std::array<std::optional<topology::Edge>, 3>> intersections(
        level_set_values.size(), std::array<std::optional<topology::Edge>, 3>() );

    const auto add_any_intersections = [&]( const topology::Edge& e, const size_t line_ii ) {
        const std::pair<double, double> bounds = {
            harmonic_func( topology::Vertex( e.dart() ) ),
            harmonic_func( topology::Vertex( phi( map, 1, e.dart() ).value() ) ) };
        for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
        {
            if( ( bounds.first < level_set_values.at( level_ii ) and
                    bounds.second >= level_set_values.at( level_ii ) ) or
                ( bounds.first >= level_set_values.at( level_ii ) and
                    bounds.second < level_set_values.at( level_ii ) ) )
            {
                intersections.at( level_ii ).at( line_ii ).emplace( e );
            }
        }
    };

    for( const topology::Edge& e : path0 ) add_any_intersections( e, 0 );
    for( const topology::Edge& e : path1 ) add_any_intersections( e, 1 );
    for( const topology::Edge& e : path2 ) add_any_intersections( e, 2 );

    std::vector<std::array<topology::Edge, 3>> out;
    out.reserve( intersections.size() - 2 );
    for( size_t level_ii = 1; level_ii < intersections.size() - 1; level_ii++ )
    {
        if( not( intersections.at( level_ii ).at( 0 ).has_value() and
                 intersections.at( level_ii ).at( 1 ).has_value() and
                 intersections.at( level_ii ).at( 2 ).has_value() ) )
            throw std::runtime_error( "Not all edges intersected" );
        out.push_back( { intersections.at( level_ii ).at( 0 ).value(),
                         intersections.at( level_ii ).at( 1 ).value(),
                         intersections.at( level_ii ).at( 2 ).value() } );
    }
    return out;
}

void nonDiskFoliations( const SweepInput& sweep,
                        const std::vector<double>& level_set_values,
                        const std::array<std::array<size_t, 3>,2>& cut_v_ids,//base then target
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

    const auto base = std::make_shared<const topology::CombinatorialMapRestriction>( bdry, keep_face_base, true );
    const topology::CombinatorialMapRestriction target( bdry, keep_face_target, true );

    const auto vertex_positions = [&sweep]( const topology::CombinatorialMap& map ) {
        const auto vertex_ids = indexingOrError( map, 0 );
        return [&sweep, vertex_ids]( const topology::Vertex& v ) -> Eigen::Vector3d {
            return sweep.mesh.points.at( vertex_ids( v ) );
        };
    };

    const std::array<std::array<topology::Vertex, 3>, 2> vs = [&]() {
        std::array<SmallVector<topology::Vertex, 3>,2> out( { SmallVector<topology::Vertex, 3>( 3, topology::Vertex() ), 
                                                              SmallVector<topology::Vertex, 3>( 3, topology::Vertex() ) } );
        const auto map_vertex_ids = indexingOrError( map, 0 );
        iterateCellsWhile( map, 0, [&]( const topology::Vertex& v ) {
            const size_t id = map_vertex_ids( v );

            auto it = std::find( cut_v_ids.at( 0 ).begin(), cut_v_ids.at( 0 ).end(), id );
            if( it != cut_v_ids.at( 0 ).end() )
            {
                out.at( 0 ).at( std::distance( cut_v_ids.at( 0 ).begin(), it ) ) = v;
            }

            it = std::find( cut_v_ids.at( 1 ).begin(), cut_v_ids.at( 1 ).end(), id );
            if( it != cut_v_ids.at( 1 ).end() )
            {
                out.at( 1 ).at( std::distance( cut_v_ids.at( 1 ).begin(), it ) ) = v;
            }
            return true;
        } );

        return std::array<std::array<topology::Vertex, 3>, 2>{ std::array<topology::Vertex, 3>{ out.at(0).at(0), out.at(0).at(1), out.at(0).at(2) },
                                                               std::array<topology::Vertex, 3>{ out.at(1).at(0), out.at(1).at(1), out.at(1).at(2) } };
    }();

    const auto harmonic_func = [vol_vertex_ids = indexingOrError( map, 0 ), &sol]( const topology::Vertex& v ) {
        return sol( vol_vertex_ids( v ) );
    };

    const std::vector<std::array<topology::Edge, 3>> tunnel_loop_intersections =
        traceIntersections( map, vs, level_set_values, vertex_positions( map ), harmonic_func );

    using namespace topology;
    std::vector<reparam::FoliationLeaf> leaves;

    { // Base level set
        const auto find_on_base = [&]( const Vertex& v ) {
            std::optional<Vertex> out;
            iterateDartsOfCell( map, v, [&]( const Dart& d ) {
                if( keep_face_base( d ) )
                {
                    out.emplace( d );
                    return false;
                }
                return true;
            } );
            return out.value();
        };
        const std::array<topology::Vertex, 3> cut_points{
            bdry.fromUnderlyingCell( find_on_base( vs.at( 0 ).at( 0 ) ) ),
            bdry.fromUnderlyingCell( find_on_base( vs.at( 0 ).at( 1 ) ) ),
            bdry.fromUnderlyingCell( find_on_base( vs.at( 0 ).at( 2 ) ) ) };
        const auto base_positions = vertex_positions( bdry );
        leaves.push_back( leafFromSphereLevelSet( base, base_positions, cut_points ) );

        std::cout << vertex_positions( map )( vs.at( 0 ).at( 0 ) ).transpose() << " vs " << base_positions( cut_points.at( 0 ) ).transpose() << std::endl;
        std::cout << vertex_positions( map )( vs.at( 0 ).at( 1 ) ).transpose() << " vs " << base_positions( cut_points.at( 1 ) ).transpose() << std::endl;
        std::cout << vertex_positions( map )( vs.at( 0 ).at( 2 ) ).transpose() << " vs " << base_positions( cut_points.at( 2 ) ).transpose() << std::endl;
    }

    for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
    { // Midway level set
        const auto level_set = std::make_shared<const topology::LevelSetCMap>( map, harmonic_func, level_set_values[level_ii] );
        const auto v_pos = vertex_positions( map );
        const auto level_set_positions = topology::levelSetVertexPositions( *level_set, v_pos );
        const auto level_set_tri = std::make_shared<const topology::DelaunayTriangulation>( level_set, level_set_positions );
        const auto tri_positions =
            topology::delaunayTriangulationVertexPositions( *level_set_tri, level_set_positions );

        const auto find_on_level = [&]( const Edge& e ) {
            std::optional<Vertex> out;
            iterateDartsOfCell( map, e, [&]( const Dart& d ) {
                if( level_set->isInMap( topology::Vertex( d ) ) )
                {
                    out.emplace( d );
                    return false;
                }
                return true;
            } );
            return out.value();
        };

        const std::array<topology::Vertex, 3> cut_points{
            find_on_level( tunnel_loop_intersections.at( level_ii - 1 ).at( 0 ).dart() ),
            find_on_level( tunnel_loop_intersections.at( level_ii - 1 ).at( 1 ).dart() ),
            find_on_level( tunnel_loop_intersections.at( level_ii - 1 ).at( 2 ).dart() ) };

        std::cout << "LEVEL " << level_ii << " at value " << level_set_values.at( level_ii ) << std::endl;
        leaves.push_back( leafFromSphereLevelSet( level_set_tri, tri_positions, cut_points ) );
    }

    { // target level set
        const auto find_on_target = [&]( const Vertex& v ) {
            std::optional<Vertex> out;
            iterateDartsOfCell( map, v, [&]( const Dart& d ) {
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
        const std::array<topology::Vertex, 3> cut_points{
            rev_map->fromUnderlyingCell( bdry.fromUnderlyingCell( find_on_target( vs.at( 1 ).at( 0 ) ) ) ),
            rev_map->fromUnderlyingCell( bdry.fromUnderlyingCell( find_on_target( vs.at( 1 ).at( 1 ) ) ) ),
            rev_map->fromUnderlyingCell( bdry.fromUnderlyingCell( find_on_target( vs.at( 1 ).at( 2 ) ) ) ) };

        leaves.push_back( leafFromSphereLevelSet( rev_map, rev_positions, cut_points ) );
    }

    callback( leaves );
}

void bulletFoliation( const std::vector<double>& level_set_values,
                      const std::function<void( const std::vector<reparam::FoliationLeaf>& )>& callback )
{
    const SweepInput sweep = SweepInputTestCases::bullet_full();

    const std::array<std::array<size_t, 3>,2> cut_v_ids({ { 2225, 98, 2237 }, { 415, 0, 27 } }); // FIXME

    nonDiskFoliations( sweep, level_set_values, cut_v_ids, callback );
}

TEST_CASE( "Level set parameterization of bullet in sphere" )
{
    const size_t n_levels = 18;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    SimplicialComplex level_sets;
    const std::vector<Eigen::Vector2d> square_points = util::generatePointsInPolygon( 150, 4 );

    std::vector<std::vector<Eigen::Vector3d>> mapped_points( square_points.size() );

    bulletFoliation( level_set_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
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
                    std::cerr << "NO VALUE level: " << level_ii << " pt: " << square_points.at( point_ii ).transpose() << std::endl;
                    break;
                }
                const auto space_pt = leaf.space_mapping->evaluate( param_pt.value().first, param_pt.value().second );
                mapped_points.at( point_ii ).push_back( space_pt );
            }

            iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
                addTriangleNoDuplicateChecking( level_sets, triangleOfFace<3>( cmap, positions, f ) );
                return true;
            } );

            const auto tutte_positions = [tp=square_mapping->vertPositions(),&level_ii]( const topology::Vertex& v ) -> Eigen::Vector3d {
                const Eigen::Vector2d tut = tp( v );
                return ( Eigen::Vector3d() << tut, double( level_ii ) ).finished();
            };

            SimplicialComplex tutte_level_set;
            iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
                addTriangleNoDuplicateChecking( tutte_level_set, triangleOfFace<3>( cmap, tutte_positions, f ) );
                return true;
            } );
            io::VTKOutputObject output( tutte_level_set );
            io::outputSimplicialFieldToVTK( output, "bullet_level_tutte_" + std::to_string( level_ii ) + ".vtu" );
        }
    } );

    io::VTKOutputObject output( level_sets );
    io::outputSimplicialFieldToVTK( output, "bullet_level_sets.vtu" );

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
    io::outputSimplicialFieldToVTK( output4, "bullet_traces.vtu" );
}