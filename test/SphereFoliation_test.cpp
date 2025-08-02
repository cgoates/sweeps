#include <catch2/catch_test_macros.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <TriangleParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <QuadMeshCombinatorialMap.hpp>
#include <TriangleMeshMapping.hpp>
#include <OrbifoldMapping.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>
#include <Foliation.hpp>
#include <LevelSetCMap.hpp>
#include <ReversedCombinatorialMap.hpp>
#include <DelaunayTriangulation.hpp>
#include <VTKOutput.hpp>
#include <MeshInput.hpp>
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
#include <MultiPatchDecomposition.hpp>
#include <HexMeshOptimization.hpp>
#include <MFEMOutput.hpp>
#include <FittingUtilities.hpp>

#include <fstream>
#include <sstream>

#include <Eigen/Dense>

constexpr bool ORBIFOLD = true;

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

    const auto cut1 = topology::shortestPath( cmap, positions, cut_vertices.at( 0 ), testEqualVertices( cmap_vert_ids, cut_vertices.at( 1 ) ) );
    topology::GlobalCellMarker cut_marker( cmap, 1 );
    for( const auto& e : cut1 ) cut_marker.mark( cmap, e );
    const auto cut2 = topology::shortestPath( cmap, [&]( const topology::Edge& e ){
        return cut_marker.isMarked( e ) ? std::numeric_limits<double>::max() : edgeLength( cmap, positions, e );
    }, cut_vertices.at( 1 ), testEqualVertices( cmap_vert_ids, cut_vertices.at( 2 ) ) );

    io::outputEdgeChain( cmap, positions, util::concatenate( cut1, cut2 ), "level_set_cut_" + std::to_string( j ) + ".vtu" );
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
        std::make_shared<Eigen::MatrixX2d>( reparam::tutteEmbedding( *cut_cmap, positions, constraints_func, reparam::Laplace2dEdgeWeights::Cotangent ) );

    const auto tutte_positions = [cut_vert_ids, tutte = *leaf.tutte]( const topology::Vertex& v ) -> Eigen::Vector2d {
        return ( Eigen::Vector2d() << tutte( cut_vert_ids( v ), 0 ), tutte( cut_vert_ids( v ), 1 ) ).finished();
    };
    const auto cut_atlas = std::make_shared<const param::TriangleParametricAtlas>( cut_cmap );
    const auto atlas = std::make_shared<const param::TriangleParametricAtlas>( cmap_ptr );
    leaf.tutte_mapping = std::make_shared<const mapping::TriangleMeshMapping>( cut_atlas, tutte_positions, 2 );
    leaf.space_mapping = std::make_shared<const mapping::TriangleMeshMapping>( atlas, positions, 3 );
    return leaf;
}

reparam::FoliationLeaf
    orbifoldLeafFromSphereLevelSet( const std::shared_ptr<const topology::CombinatorialMap>& cmap_ptr,
                                    const VertexPositionsFunc& positions,
                                    const std::array<topology::Vertex, 3>& cut_vertices )
{
    const auto& cmap = *cmap_ptr;
    const auto cmap_vert_ids = indexingOrError( cmap, 0 );
    std::set<topology::Cell> cuts;

    const auto cut1 = topology::shortestPath(
        cmap, positions, cut_vertices.at( 0 ), testEqualVertices( cmap_vert_ids, cut_vertices.at( 1 ) ) );
    topology::GlobalCellMarker cut_marker( cmap, 1 );
    for( const auto& e : cut1 ) cut_marker.mark( cmap, e );
    const auto cut2 = topology::shortestPath(
        cmap,
        [&]( const topology::Edge& e ) {
            return cut_marker.isMarked( e ) ? std::numeric_limits<double>::max() : edgeLength( cmap, positions, e );
        },
        cut_vertices.at( 1 ),
        testEqualVertices( cmap_vert_ids, cut_vertices.at( 2 ) ) );

    io::outputEdgeChain(
        cmap, positions, util::concatenate( cut1, cut2 ), "level_set_cut_" + std::to_string( j ) + ".vtu" );
    j++;
    cuts.insert( cut1.begin(), cut1.end() );
    cuts.insert( cut2.begin(), cut2.end() );

    const auto cut_cmap = std::make_shared<const topology::CutCombinatorialMap>( cmap, cuts );
    const auto cut_vert_ids = indexingOrError( *cut_cmap, 0 );

    reparam::FoliationLeaf leaf{};
    leaf.tutte =
        std::make_shared<Eigen::MatrixX2d>( reparam::tutteOrbifoldEmbedding( *cut_cmap, positions, cut_vertices, reparam::Laplace2dEdgeWeights::Cotangent ) );

    const auto tutte_positions = [cut_vert_ids, tutte = *leaf.tutte]( const topology::Vertex& v ) -> Eigen::Vector2d {
        return ( Eigen::Vector2d() << tutte( cut_vert_ids( v ), 0 ), tutte( cut_vert_ids( v ), 1 ) ).finished();
    };
    const auto cut_atlas = std::make_shared<const param::TriangleParametricAtlas>( cut_cmap );
    const auto atlas = std::make_shared<const param::TriangleParametricAtlas>( cmap_ptr );
    leaf.tutte_mapping = std::make_shared<const mapping::OrbifoldMapping>( cut_atlas, tutte_positions );
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
        io::outputEdgeChain( map, map_positions, path0, "sphere_cut_lines0.vtu" );
        io::outputEdgeChain( map, map_positions, path1, "sphere_cut_lines1.vtu" );
        io::outputEdgeChain( map, map_positions, path2, "sphere_cut_lines2.vtu" );
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
            if( not out.has_value() )
                FAIL( "Vertex " << v << " at position " << vertex_positions( map )( v ).transpose()
                                << " not found on base" );
            return out.value();
        };
        const std::array<topology::Vertex, 3> cut_points{
            bdry.fromUnderlyingCell( find_on_base( vs.at( 0 ).at( 0 ) ) ),
            bdry.fromUnderlyingCell( find_on_base( vs.at( 0 ).at( 1 ) ) ),
            bdry.fromUnderlyingCell( find_on_base( vs.at( 0 ).at( 2 ) ) ) };
        const auto base_positions = vertex_positions( bdry );
        if( ORBIFOLD )
            leaves.push_back( orbifoldLeafFromSphereLevelSet( base, base_positions, cut_points ) );
        else
            leaves.push_back( leafFromSphereLevelSet( base, base_positions, cut_points ) );
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
        if( ORBIFOLD )
            leaves.push_back( orbifoldLeafFromSphereLevelSet( level_set_tri, tri_positions, cut_points ) );
        else
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

        if( ORBIFOLD )
            leaves.push_back( orbifoldLeafFromSphereLevelSet( rev_map, rev_positions, cut_points ) );
        else
            leaves.push_back( leafFromSphereLevelSet( rev_map, rev_positions, cut_points ) );
    }

    callback( leaves );
}

void bulletFoliation( const std::vector<double>& level_set_values,
                      const std::function<void( const std::vector<reparam::FoliationLeaf>& )>& callback )
{
    const SweepInput sweep = SweepInputTestCases::bullet_full();

    const std::array<std::array<size_t, 3>,2> cut_v_ids = {{ { 2225, 98, 2237 }, { 415, 0, 27 } }};

    nonDiskFoliations( sweep, level_set_values, cut_v_ids, callback );
}

void capsuleFoliation( const std::vector<double>& level_set_values,
                      const std::function<void( const std::vector<reparam::FoliationLeaf>& )>& callback )
{
    const SweepInput sweep = SweepInputTestCases::capsule();

    const std::array<std::array<size_t, 3>,2> cut_v_ids = {{ { 4989, 149, 5475 }, { 1, 785, 0 } }};

    nonDiskFoliations( sweep, level_set_values, cut_v_ids, callback );
}

TEST_CASE( "Level set parameterization of bullet in sphere" )
{
    const size_t n_levels = 18;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    SimplicialComplex level_sets;

    const std::vector<Eigen::Vector2d> square_points = [&]() {
        if( not ORBIFOLD )
            return util::generatePointsInPolygon( 150, 4 );

        std::vector<Eigen::Vector2d> square_points;
        for( size_t i = 1; i < 10; i++ )
        {
            for( size_t j = 0; j < 10; j++ )
            {
                square_points.push_back( Eigen::Vector2d( double( i ) / 10.0, double( j ) / 10.0 ) );
            }
        }
        return square_points;
    }();

    std::vector<std::vector<Eigen::Vector3d>> mapped_points( square_points.size() );

    const std::string output_prefix = ORBIFOLD ? "bulleto" : "bullet";
    bulletFoliation( level_set_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
        std::cout << "-------------------------" << std::endl;
        for( size_t level_ii = 0; level_ii < leaves.size(); level_ii++ )
        {
            std::cout << "LEVEL " << level_ii << " at value " << level_set_values.at( level_ii ) << std::endl;
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
            if( true )
            {
                SimplicialComplex individual_level_set;
                addAllTriangles( individual_level_set, leaf.tutte_mapping->parametricAtlas().cmap(),
                                                                            leaf.space_mapping->vertPositions() );

                io::VTKOutputObject output( individual_level_set );
                Eigen::MatrixX2d tutte( individual_level_set.points.size(), 2 );
                Eigen::Index row = 0;
                iterateCellsWhile( leaf.tutte_mapping->parametricAtlas().cmap(), 0, [&]( const auto& vert ) {
                    tutte.row( row++ ) = leaf.tutte_mapping->vertPositions()( vert );
                    return true;
                } );

                output.addVertexField( "tutte", tutte );
                io::outputSimplicialFieldToVTK( output, output_prefix + "_level_set_" + std::to_string( level_ii ) + ".vtu" );
            }

            const auto tutte_positions = [tp=square_mapping->vertPositions()]( const topology::Vertex& v ) -> Eigen::Vector3d {
                const Eigen::Vector2d tut = tp( v );
                return ( Eigen::Vector3d() << tut, 0.0 ).finished();
            };

            SimplicialComplex tutte_level_set;
            iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
                addTriangleNoDuplicateChecking( tutte_level_set, triangleOfFace<3>( cmap, tutte_positions, f ) );
                return true;
            } );
            io::VTKOutputObject output( tutte_level_set );
            io::outputSimplicialFieldToVTK( output, output_prefix + "_level_tutte_" + std::to_string( level_ii ) + ".vtu" );

            if( ORBIFOLD )
            {
                SimplicialComplex tutte_level_set2;
                const auto& square_mapping2 = static_cast<const mapping::OrbifoldMapping&>( *leaf.tutte_mapping );
                square_mapping2.iterateTriangles( [&]( const Triangle<2>& t ) {
                    const Triangle<3> t3{ (Eigen::Vector3d() << t.v1, 0 ).finished(),
                                           (Eigen::Vector3d() << t.v2, 0 ).finished(),
                                           (Eigen::Vector3d() << t.v3, 0 ).finished() };
                    addTriangleNoDuplicateChecking( tutte_level_set2, t3 );
                    return true;
                } );

                io::VTKOutputObject output2( tutte_level_set2 );
                io::outputSimplicialFieldToVTK( output2, output_prefix + "_level_tutte2_" + std::to_string( level_ii ) + ".vtu" );
            }
        }
    } );

    io::VTKOutputObject output( level_sets );
    io::outputSimplicialFieldToVTK( output, output_prefix + "_level_sets.vtu" );

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
    io::outputSimplicialFieldToVTK( output4, output_prefix + "_traces.vtu" );
}

double checkJacobians( const basis::SplineSpace& ss, const Eigen::MatrixXd& cpts )
{
    double min_val = std::numeric_limits<double>::infinity();
    double max_val = -std::numeric_limits<double>::infinity();

    eval::SplineSpaceEvaluator evaler( ss, 1 );

    min_val = std::numeric_limits<double>::infinity();
    max_val = -std::numeric_limits<double>::infinity();

    iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const topology::Volume& vol ) {
        evaler.localizeElement( vol );
        iterateAdjacentCellsOfRestrictedCell( ss.basisComplex().parametricAtlas().cmap(), vol, 2, 0, [&]( const topology::Vertex& v ) {
            evaler.localizePoint( ss.basisComplex().parametricAtlas().parentPoint( v ) );
            const double val = evaler.evaluateJacobian( cpts.transpose() ).determinant();
            min_val = std::min( min_val, val );
            max_val = std::max( max_val, val );
            return true;
        } );
        return true;
    } );

    std::cout << "Min jacobian: " << min_val << std::endl;
    std::cout << "Max jacobian: " << max_val << std::endl;
    return min_val;
}

TEST_CASE( "Hex mesh on the bullet in sphere" )
{
    const size_t n_levels = 18;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    // Load the quad mesh
    const auto pr = io::loadOBJFile( SRC_HOME "/test/data/bullet.obj" );
    const auto& faces = pr.first;
    const auto& points = pr.second;
    const topology::QuadMeshCombinatorialMap quad_layout = [&](){
        std::vector<std::array<VertexId, 4>> quads;
        quads.reserve( faces.size() );
        for( const auto& f : faces )
        {
            if( f.size() == 4 )
            {
                quads.push_back( { VertexId( f.at( 0 ) ), VertexId( f.at( 1 ) ), VertexId( f.at( 2 ) ), VertexId( f.at( 3 ) ) } );
            }
            else
            {
                throw std::runtime_error( "Only quads are supported in this test." );
            }
        }
        return topology::QuadMeshCombinatorialMap( quads, points.size() );
    }();

    const topology::MultiPatchDecomposition decomp = topology::multiPatchDecomposition( quad_layout );
    const topology::MultiPatchCombinatorialMap cmap2d( decomp.constituents, decomp.connections );

    topology::GlobalCellMarker m( quad_layout, 0 );
    std::map<topology::Dart::IndexType, Eigen::Vector3d> cmap2d_positions;

    const auto quad_layout_vids = indexingOrError( quad_layout, 0 );
    synchronizedFlood2d( quad_layout,
                         cmap2d,
                         decomp.unstructured_first_corner,
                         topology::Dart( 0 ),
                         [&]( const topology::Dart& d1, const topology::Dart& d2 ) {
                             if( not m.isMarked( topology::Vertex( d1 ) ) )
                             {
                                 m.mark( quad_layout, topology::Vertex( d1 ) );
                                 const auto vid = quad_layout_vids( topology::Vertex( d1 ) );
                                 cmap2d_positions.insert(
                                     { lowestDartId( cmap2d, topology::Vertex( d2 ) ), points.at( vid ) } );
                             }
                             return true;
                         } );

    // Create the 3d linear spline space
    const basis::MultiPatchSplineSpace ss = [&]() {
        std::vector<std::shared_ptr<const basis::TPSplineSpace>> constituents;
        constituents.reserve( cmap2d.constituents().size() );

        const basis::KnotVector kv_u = basis::unitIntervalKnotVectorWithNElems( level_set_values.size() - 1, 1 );

        for( const auto& constituent : cmap2d.constituents() )
        {
            const auto& components = tensorProductComponentCMaps( *constituent );
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( cellCount( *components.at( 0 ), 1 ), 1 );
            const basis::KnotVector kv_t = basis::unitIntervalKnotVectorWithNElems( cellCount( *components.at( 1 ), 1 ), 1 );

            constituents.push_back( std::make_shared<const basis::TPSplineSpace>( basis::buildBSpline( {kv_s, kv_t, kv_u }, {1, 1, 1} ) ) );
        }
        return buildMultiPatchSplineSpace( constituents, topology::connectionsOfSweptMultipatch( cmap2d.connections() ) );
    }();

    const auto sweep_to_source_vertex = [&]( const topology::Vertex& v ) {
        const auto [patch_id, d_patch] = ss.basisComplex().parametricAtlas().cmap().toLocalDart( v.dart() );
        const auto unflattened_cell = unflattenCell( ss.subSpaces().at( patch_id )->basisComplex().parametricAtlas().cmap(), topology::Vertex( d_patch ) );
        const topology::Vertex v2d( lowestDartId(
            cmap2d, topology::Vertex( cmap2d.toGlobalDart( patch_id, unflattened_cell.first.value().dart() ) ) ) );
        const size_t leaf_ii = unflattened_cell.second.has_value() ? unflattened_cell.second.value().dart().id() : level_set_values.size() - 1;

        return std::pair<topology::Vertex, size_t>( v2d, leaf_ii );
    };

    const auto get_tutte_points = [&]( const reparam::FoliationLeaf& leaf ) {
        std::map<topology::Dart::IndexType, Eigen::Vector2d> tutte_points;
        // Find the tutte domain coordinates of the source cmap vertices.

        for( const auto& pr : cmap2d_positions )
        {
            const auto [cell, parent_pt] = leaf.space_mapping->closestPoint( pr.second );
            const Eigen::Vector2d tutte_pt = leaf.tutte_mapping->evaluate( cell, parent_pt );
            tutte_points.insert( { pr.first, tutte_pt } );
        }
        
        return tutte_points;
    };

    const SweepInput sweep = SweepInputTestCases::bullet_full();

    const std::array<std::array<size_t, 3>,2> cut_v_ids = {{ { 2225, 98, 2237 }, { 415, 0, 27 } }};

    nonDiskFoliations( sweep, level_set_values, cut_v_ids, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
        const auto tutte_points = get_tutte_points( leaves.front() );
        const Eigen::MatrixXd fit_cpts = reparam::fitLinearMeshToLeaves( ss, leaves, [&]( const topology::Vertex& v ){
            const auto [v2d, leaf_ii] = sweep_to_source_vertex( v );
            return std::pair<Eigen::Vector2d, size_t>( tutte_points.at( v2d.dart().id() ), leaf_ii );
        } );

        io::outputBezierMeshToVTK( ss, fit_cpts, "fit_to_bullet_multi_patch.vtu" );

        checkJacobians( ss, fit_cpts );

        // Optimize the mesh
        {
            const topology::TetMeshCombinatorialMap tet_map( sweep.mesh );
            const auto tet_vert_idx = indexingOrError( tet_map, 0 );
            const auto tet_vertex_positions = [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
                return sweep.mesh.points.at( tet_vert_idx( v ) );
            };
            const auto hex_vertex_positions = eval::vertexPositionsFromManifold( ss, fit_cpts.transpose() );
            const auto new_positions = fitting::optimizeMesh( tet_map, tet_vertex_positions, ss.basisComplex().parametricAtlas().cmap(), hex_vertex_positions,
                                                    "output.vtk" );

            const Eigen::MatrixXd new_fit_cpts = fitting::linearControlPointsFromVertexPositions( ss, new_positions );

            CHECK( checkJacobians( ss, new_fit_cpts ) > 0.0 );
            io::outputBezierMeshToVTK( ss, new_fit_cpts, "fit_to_bullet_multi_patch_optimized.vtu" );
        }
        
    } );
}

TEST_CASE( "Level set parameterization of capsule in ellipsoid" )
{
    const size_t n_levels = 18;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    SimplicialComplex level_sets;

    const std::vector<Eigen::Vector2d> square_points = [&]() {
        if( not ORBIFOLD )
            return util::generatePointsInPolygon( 150, 4 );

        std::vector<Eigen::Vector2d> square_points;
        for( size_t i = 1; i < 10; i++ )
        {
            for( size_t j = 0; j < 10; j++ )
            {
                square_points.push_back( Eigen::Vector2d( double( i ) / 10.0, double( j ) / 10.0 ) );
            }
        }
        return square_points;
    }();

    std::vector<std::vector<Eigen::Vector3d>> mapped_points( square_points.size() );

    const std::string output_prefix = ORBIFOLD ? "capsule" : "capsule_noorbifold";
    capsuleFoliation( level_set_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
        std::cout << "-------------------------" << std::endl;
        for( size_t level_ii = 0; level_ii < leaves.size(); level_ii++ )
        {
            std::cout << "LEVEL " << level_ii << " at value " << level_set_values.at( level_ii ) << std::endl;
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
            if( true )
            {
                SimplicialComplex individual_level_set;
                addAllTriangles( individual_level_set, leaf.tutte_mapping->parametricAtlas().cmap(),
                                                                            leaf.space_mapping->vertPositions() );

                io::VTKOutputObject output( individual_level_set );
                Eigen::MatrixX2d tutte( individual_level_set.points.size(), 2 );
                Eigen::Index row = 0;
                iterateCellsWhile( leaf.tutte_mapping->parametricAtlas().cmap(), 0, [&]( const auto& vert ) {
                    tutte.row( row++ ) = leaf.tutte_mapping->vertPositions()( vert );
                    return true;
                } );

                output.addVertexField( "tutte", tutte );
                io::outputSimplicialFieldToVTK( output, output_prefix + "_level_set_" + std::to_string( level_ii ) + ".vtu" );
            }

            const auto tutte_positions = [tp=square_mapping->vertPositions()]( const topology::Vertex& v ) -> Eigen::Vector3d {
                const Eigen::Vector2d tut = tp( v );
                return ( Eigen::Vector3d() << tut, 0.0 ).finished();
            };

            SimplicialComplex tutte_level_set;
            iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
                addTriangleNoDuplicateChecking( tutte_level_set, triangleOfFace<3>( cmap, tutte_positions, f ) );
                return true;
            } );
            io::VTKOutputObject output( tutte_level_set );
            io::outputSimplicialFieldToVTK( output, output_prefix + "_level_tutte_" + std::to_string( level_ii ) + ".vtu" );

            if( ORBIFOLD )
            {
                SimplicialComplex tutte_level_set2;
                const auto& square_mapping2 = static_cast<const mapping::OrbifoldMapping&>( *leaf.tutte_mapping );
                square_mapping2.iterateTriangles( [&]( const Triangle<2>& t ) {
                    const Triangle<3> t3{ (Eigen::Vector3d() << t.v1, 0 ).finished(),
                                           (Eigen::Vector3d() << t.v2, 0 ).finished(),
                                           (Eigen::Vector3d() << t.v3, 0 ).finished() };
                    addTriangleNoDuplicateChecking( tutte_level_set2, t3 );
                    return true;
                } );

                io::VTKOutputObject output2( tutte_level_set2 );
                io::outputSimplicialFieldToVTK( output2, output_prefix + "_level_tutte2_" + std::to_string( level_ii ) + ".vtu" );
            }
        }
    } );

    io::VTKOutputObject output( level_sets );
    io::outputSimplicialFieldToVTK( output, output_prefix + "_level_sets.vtu" );

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
    io::outputSimplicialFieldToVTK( output4, output_prefix + "_traces.vtu" );
}

TEST_CASE( "Hex mesh on the capsule" )
{
    constexpr bool LOAD_FROM_FILE = true;
    const size_t n_levels = 10;
    const std::vector<double> level_set_values = util::linspace( 0, 1, n_levels );

    // Load the quad mesh
    const auto pr = io::loadOBJFile( SRC_HOME "/test/data/capsule.obj" );
    const auto& faces = pr.first;
    const auto& points = pr.second;
    const topology::QuadMeshCombinatorialMap quad_layout = [&](){
        std::vector<std::array<VertexId, 4>> quads;
        quads.reserve( faces.size() );
        for( const auto& f : faces )
        {
            if( f.size() == 4 )
            {
                quads.push_back( { VertexId( f.at( 0 ) ), VertexId( f.at( 1 ) ), VertexId( f.at( 2 ) ), VertexId( f.at( 3 ) ) } );
            }
            else
            {
                throw std::runtime_error( "Only quads are supported in this test." );
            }
        }
        return topology::QuadMeshCombinatorialMap( quads, points.size() );
    }();

    const topology::MultiPatchDecomposition decomp = topology::multiPatchDecomposition( quad_layout );
    const topology::MultiPatchCombinatorialMap cmap2d( decomp.constituents, decomp.connections );

    topology::GlobalCellMarker m( quad_layout, 0 );
    std::map<topology::Dart::IndexType, Eigen::Vector3d> cmap2d_positions;

    const auto quad_layout_vids = indexingOrError( quad_layout, 0 );
    synchronizedFlood2d( quad_layout,
                         cmap2d,
                         decomp.unstructured_first_corner,
                         topology::Dart( 0 ),
                         [&]( const topology::Dart& d1, const topology::Dart& d2 ) {
                             if( not m.isMarked( topology::Vertex( d1 ) ) )
                             {
                                 m.mark( quad_layout, topology::Vertex( d1 ) );
                                 const auto vid = quad_layout_vids( topology::Vertex( d1 ) );
                                 cmap2d_positions.insert(
                                     { lowestDartId( cmap2d, topology::Vertex( d2 ) ), points.at( vid ) } );
                             }
                             return true;
                         } );

    // Create the 3d linear spline space
    const basis::MultiPatchSplineSpace ss = [&]() {
        std::vector<std::shared_ptr<const basis::TPSplineSpace>> constituents;
        constituents.reserve( cmap2d.constituents().size() );

        const basis::KnotVector kv_u = basis::unitIntervalKnotVectorWithNElems( level_set_values.size() - 1, 1 );

        for( const auto& constituent : cmap2d.constituents() )
        {
            const auto& components = tensorProductComponentCMaps( *constituent );
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( cellCount( *components.at( 0 ), 1 ), 1 );
            const basis::KnotVector kv_t = basis::unitIntervalKnotVectorWithNElems( cellCount( *components.at( 1 ), 1 ), 1 );

            constituents.push_back( std::make_shared<const basis::TPSplineSpace>( basis::buildBSpline( {kv_s, kv_t, kv_u }, {1, 1, 1} ) ) );
        }
        return buildMultiPatchSplineSpace( constituents, topology::connectionsOfSweptMultipatch( cmap2d.connections() ) );
    }();

    Eigen::MatrixXd final_fit;

    if( LOAD_FROM_FILE )
    {
        final_fit = readFromFile( SRC_HOME "/test/data/capsule_control_points.txt" );
    }
    else
    {
        const auto sweep_to_source_vertex = [&]( const topology::Vertex& v ) {
            const auto [patch_id, d_patch] = ss.basisComplex().parametricAtlas().cmap().toLocalDart( v.dart() );
            const auto unflattened_cell = unflattenCell( ss.subSpaces().at( patch_id )->basisComplex().parametricAtlas().cmap(), topology::Vertex( d_patch ) );
            const topology::Vertex v2d( lowestDartId(
                cmap2d, topology::Vertex( cmap2d.toGlobalDart( patch_id, unflattened_cell.first.value().dart() ) ) ) );
            const size_t leaf_ii = unflattened_cell.second.has_value() ? unflattened_cell.second.value().dart().id() : level_set_values.size() - 1;

            return std::pair<topology::Vertex, size_t>( v2d, leaf_ii );
        };

        const auto get_tutte_points = [&]( const reparam::FoliationLeaf& leaf ) {
            std::map<topology::Dart::IndexType, Eigen::Vector2d> tutte_points;
            // Find the tutte domain coordinates of the source cmap vertices.

            for( const auto& pr : cmap2d_positions )
            {
                const auto [cell, parent_pt] = leaf.space_mapping->closestPoint( pr.second );
                const Eigen::Vector2d tutte_pt = leaf.tutte_mapping->evaluate( cell, parent_pt );
                tutte_points.insert( { pr.first, tutte_pt } );
            }
            
            return tutte_points;
        };

        const SweepInput sweep = SweepInputTestCases::capsule();

        const std::array<std::array<size_t, 3>,2> cut_v_ids = {{ { 4989, 149, 5475 }, { 1, 785, 0 } }};

        nonDiskFoliations( sweep, level_set_values, cut_v_ids, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
            const auto tutte_points = get_tutte_points( leaves.front() );
            const Eigen::MatrixXd fit_cpts = reparam::fitLinearMeshToLeaves( ss, leaves, [&]( const topology::Vertex& v ){
                const auto [v2d, leaf_ii] = sweep_to_source_vertex( v );
                return std::pair<Eigen::Vector2d, size_t>( tutte_points.at( v2d.dart().id() ), leaf_ii );
            } );

            io::outputBezierMeshToVTK( ss, fit_cpts, "fit_to_capsule_multi_patch.vtu" );

            double min_val = std::numeric_limits<double>::infinity();
            double max_val = -std::numeric_limits<double>::infinity();

            SimplicialComplex inverted_points;
            std::vector<double> jacobian_values;

            eval::SplineSpaceEvaluator evaler( ss, 1 );
            iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const topology::Volume& vol ) {
                evaler.localizeElement( vol );
                iterateAdjacentCellsOfRestrictedCell( ss.basisComplex().parametricAtlas().cmap(), vol, 2, 0, [&]( const topology::Vertex& v ) {
                    evaler.localizePoint( ss.basisComplex().parametricAtlas().parentPoint( v ) );
                    const double val = evaler.evaluateJacobian( fit_cpts.transpose() ).determinant();
                    min_val = std::min( min_val, val );
                    max_val = std::max( max_val, val );
                    if( val < 0 )
                    {
                        inverted_points.simplices.emplace_back( inverted_points.points.size() );
                        inverted_points.points.push_back( evaler.evaluateManifold( fit_cpts.transpose() ) );
                        jacobian_values.push_back( val );
                    }
                    return true;
                } );
                return true;
            } );

            std::cout << "Min jacobian: " << min_val << std::endl;
            std::cout << "Max jacobian: " << max_val << std::endl;

            io::VTKOutputObject inverted_points_output( inverted_points );
            inverted_points_output.addVertexField( "jacobian", Eigen::Map<Eigen::VectorXd>( jacobian_values.data(), jacobian_values.size() ) );
            io::outputSimplicialFieldToVTK( inverted_points_output, "inverted_points_capsule.vtu" );
            
            // Optimize the mesh
            {
                const topology::TetMeshCombinatorialMap tet_map( sweep.mesh );
                const auto tet_vert_idx = indexingOrError( tet_map, 0 );
                const auto tet_vertex_positions = [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
                    return sweep.mesh.points.at( tet_vert_idx( v ) );
                };
                const auto hex_vertex_positions = eval::vertexPositionsFromManifold( ss, fit_cpts.transpose() );
                const auto new_positions = fitting::optimizeMesh( tet_map, tet_vertex_positions, ss.basisComplex().parametricAtlas().cmap(), hex_vertex_positions,
                                                        "output.vtk" );

                const Eigen::MatrixXd new_fit_cpts = fitting::linearControlPointsFromVertexPositions( ss, new_positions );

                final_fit = new_fit_cpts;

                writeToFile( final_fit, SRC_HOME "/test/data/capsule_control_points.txt" );
            }
        } );
    }

    checkJacobians( ss, final_fit );

    io::outputBezierMeshToVTK( ss, final_fit, "fit_to_capsule_multi_patch_optim.vtu" );

    const basis::MultiPatchSplineSpace quadratic_ss = degreeRefineOrCoarsen(
        ss,
        [&, kv_u = basis::unitIntervalKnotVectorWithNElems( level_set_values.size() - 1, 2 )]( const size_t patch_ii ) {
            const auto& components = tensorProductComponentCMaps( *cmap2d.constituents().at( patch_ii ) );
            const basis::KnotVector kv_s =
                basis::unitIntervalKnotVectorWithNElems( cellCount( *components.at( 0 ), 1 ), 2 );
            const basis::KnotVector kv_t =
                basis::unitIntervalKnotVectorWithNElems( cellCount( *components.at( 1 ), 1 ), 2 );
            return basis::DegreeAndKnotVector( { { 2, 2, 2 }, { kv_s, kv_t, kv_u } } );
        } );

    const Eigen::MatrixXd quadratic_fit_cpts = fitting::fitToManifold( ss, final_fit, quadratic_ss );

    checkJacobians( quadratic_ss, quadratic_fit_cpts );

    io::outputBezierMeshToVTK( quadratic_ss, quadratic_fit_cpts, "fit_to_capsule_multi_patch_quadratic.vtu" );

    io::outputMultiPatchSplinesToMFEM( quadratic_ss, quadratic_fit_cpts.transpose(), "optimized_fit_to_capsule_multi_patch.mesh" );
}