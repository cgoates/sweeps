#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>
#include <Logging.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <TriMeshCombinatorialMap.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CommonUtils.hpp>
#include <CutCombinatorialMap.hpp>
#include <CustomCombinatorialMap.hpp>
#include <VTKOutput.hpp>
#include <Dijkstra.hpp>
#include <CombinatorialMapRestriction.hpp>

TEST_CASE( "Laplace patch test", "" )
{
    SweepInput sweep_input = SweepInputTestCases::twelveTetCube();
    sweep_input.mesh.points.back() = Eigen::Vector3d( 0.37, 0.49, 0.55 );

    topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = reparam::sweepEmbedding(
        map, sweep_input.zero_bcs, sweep_input.one_bcs, normals, reparam::Laplace3dEdgeWeights::Cotangent );

    REQUIRE( util::equals( sol( 8 ), sweep_input.mesh.points.back()( 2 ), sweep_input.mesh.points.back()( 2 ) * 1e-15 ) );
}

TEST_CASE( "Tutte embedding patch test" )
{
    // A perterbation of this:
    // *-------*
    // | \   / |
    // |  \ /  |
    // |   *   |
    // |  / \  |
    // | /   \ |
    // *-------*
    SimplicialComplex mesh{
        { { 0, 1, 4 }, { 0, 4, 2 }, { 4, 1, 3 }, { 4, 3, 2 } },
        { { 0.0, 0.0, 2.0 }, { 1.0, 0.0, 1.2 }, { 0.0, 1.0, 0.7 }, { 1.5, 1.2, 0.2 }, { 0.7, 0.6, 0.0 } } };

    const topology::TriMeshCombinatorialMap map( mesh );

    std::vector<Eigen::Vector2d> constraints{ { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { 1.0, 1.0 } };

    const auto vids = indexingOrError( map, 0 );

    const auto constraints_func = [&]( const topology::Vertex& v ) -> std::optional<Eigen::Vector2d> {
        if( not boundaryAdjacent( map, v ) ) return std::nullopt;
        return constraints.at( vids( v ) );
    };

    const auto vert_positions = [&]( const topology::Vertex& v ){
        return mesh.points.at( vids( v ) );
    };

    const Eigen::MatrixX2d tutte = reparam::tutteEmbedding( map, vert_positions, constraints_func, reparam::Laplace2dEdgeWeights::Uniform );

    const Eigen::MatrixX2d expected = ( Eigen::MatrixX2d( 5, 2 ) << 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 0.5 ).finished();


    std::cout << tutte << std::endl;
    REQUIRE( tutte.rows() == 5 );
    for( Eigen::Index i = 0; i < 5; i++ )
    {
        CHECK( tutte( i, 0 ) == expected( i, 0 ) );
        CHECK( tutte( i, 1 ) == expected( i, 1 ) );
    }
}

TEST_CASE( "Tutte Orbifold embedding" )
{
    //              * 2
    //             /|\.
    //            / | \.
    //           /  |  \.
    //          /   |   \.
    //         /    |    \.
    //        /     |     * 1
    //       /      |    / \.
    //      /       |   /   \.
    //     /        |  /     \.
    //    /         | /       \.
    //   /          |/         \.
    //  *-----------*-----------*
    // 4            0            3
    SimplicialComplex mesh{
        { { 0, 3, 1 }, { 0, 1, 3 }, { 0, 1, 2 }, { 0, 2, 1 }, { 0, 2, 4 }, { 0, 4, 2 } },
        { { 0, 0, 0 }, {0.5, 0.5, 0}, { 0, 1, 0 }, {1, 0, 0 }, { -1, 0, 0 } }
    };
    using namespace topology;
    const CustomCombinatorialMap map( 18,
                                      2,
                                      { { 1, 2, 0, 4, 5, 3, 7, 8, 6, 10, 11, 9, 13, 14, 12, 16, 17, 15 },
                                        { 11, 10, 3, 2, 13, 6, 5, 16, 15, 14, 1, 0, 17, 4, 9, 8, 7, 12 } },
                                      { { 0, { 0, 3, 1, 0, 1, 2, 0, 2, 4, 0, 1, 3, 0, 2, 1, 0, 4, 2 } } } );
    const auto indexing = indexingOrError( map, 0 );
    const topology::CutCombinatorialMap cut_map( map, { Edge( 1 ), Edge( 4 ), Edge( 7 ) } );
    const auto cut_indexing = indexingOrError( cut_map, 0 );

    const Eigen::MatrixX2d tutte = reparam::tutteOrbifoldEmbedding( cut_map, [&]( const topology::Vertex& v ) {
        return mesh.points.at( indexing( v ) );
    }, { Vertex( 1 ), Vertex( 7 ), Vertex( 8 ) }, reparam::Laplace2dEdgeWeights::Uniform );

    std::cout << tutte.row( 3 ) << std::endl;

    io::outputCMap( cut_map, [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
        return Eigen::Vector3d( tutte( cut_indexing( v ), 0 ), tutte( cut_indexing( v ), 1 ), 0 );
    }, "tutte_orbi.vtu" );
}

std::function<bool( const topology::Vertex& )> testEqualVertices( const topology::IndexingFunc& vert_ids,
                                                                  const topology::Vertex& end_v )
{
    return [&]( const topology::Vertex& test_v ) { return vert_ids( test_v ) == vert_ids( end_v ); };
}

TEST_CASE( "Large sphere tutte orbifold embedding" )
{
    const SweepInput sweep = SweepInputTestCases::bullet_full();
    const std::array<size_t, 3> cut_vert_ids( { 415, 0, 27 } );

    const topology::TetMeshCombinatorialMap map( sweep.mesh );
    const auto indexing = indexingOrError( map, 0 );

    const topology::CombinatorialMapBoundary bdry( map );
    const auto bdry_vert_ids = indexingOrError( bdry, 0 );
    const auto keep_face_target = [&]( const topology::Face& f ) {
        return sweep.one_bcs.at( bdry_vert_ids( topology::Vertex( f.dart() ) ) ) and
               sweep.one_bcs.at( bdry_vert_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
               sweep.one_bcs.at( bdry_vert_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
    };
    const topology::CombinatorialMapRestriction target( bdry, keep_face_target, true );

    io::outputCMap( target, [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
        return sweep.mesh.points.at( bdry_vert_ids( v ) );
    }, "bdry.vtu" );


    std::array<topology::Vertex, 3> cut_vertices;
    iterateCellsWhile( target, 0, [&]( const topology::Vertex& v ) {
        if( bdry_vert_ids( v ) == cut_vert_ids.at( 0 ) ) cut_vertices.at( 0 ) = v;
        else if( bdry_vert_ids( v ) == cut_vert_ids.at( 1 ) ) cut_vertices.at( 1 ) = v;
        else if( bdry_vert_ids( v ) == cut_vert_ids.at( 2 ) ) cut_vertices.at( 2 ) = v;
        return true;
    } );

    const auto positions = [&]( const topology::Vertex& v ) {
        return sweep.mesh.points.at( bdry_vert_ids( v ) );
    };

    const auto cut1 = topology::shortestPath( target, positions, cut_vertices.at( 0 ), testEqualVertices( bdry_vert_ids, cut_vertices.at( 1 ) ) );
    topology::GlobalCellMarker cut_marker( bdry, 1 );
    for( const auto& e : cut1 ) cut_marker.mark( bdry, e );
    const auto cut2 = topology::shortestPath( target, [&]( const topology::Edge& e ){
        return cut_marker.isMarked( e ) ? std::numeric_limits<double>::max() : edgeLength( bdry, positions, e );
    }, cut_vertices.at( 1 ), testEqualVertices( bdry_vert_ids, cut_vertices.at( 2 ) ) );

    io::outputEdges( bdry, positions, util::concatenate( cut1, cut2 ), "level_set_cut.vtu" );

    std::set<topology::Cell> cuts;
    cuts.insert( cut1.begin(), cut1.end() );
    cuts.insert( cut2.begin(), cut2.end() );
    const topology::CutCombinatorialMap cut_cmap( target, cuts );
    const auto cutmap_vert_ids = indexingOrError( cut_cmap, 0 );

    const Eigen::MatrixX2d tutte = reparam::tutteOrbifoldEmbedding( cut_cmap, [&]( const topology::Vertex& v ) {
        return sweep.mesh.points.at( bdry_vert_ids( v ) );
    }, { cut_vertices.at( 0 ), cut_vertices.at( 1 ), cut_vertices.at( 2 ) }, reparam::Laplace2dEdgeWeights::InverseLength );

    io::outputCMap( cut_cmap, [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
        return Eigen::Vector3d( tutte( cutmap_vert_ids( v ), 0 ), tutte( cutmap_vert_ids( v ), 1 ), 0 );
    }, "tutte_orbi_large.vtu" );
}


TEST_CASE( "Sphere tutte orbifold embedding" )
{
    const SweepInput sweep = io::loadINPFile( SRC_HOME "/test/data/Sphere.inp", "lala", "lala" );
    const std::array<size_t, 3> cut_vert_ids( { 0, 7, 1 } );

    const topology::TetMeshCombinatorialMap map( sweep.mesh );
    const auto indexing = indexingOrError( map, 0 );

    const topology::CombinatorialMapBoundary bdry( map );
    const topology::CombinatorialMapRestriction target( bdry, [&]( const auto& ) { return true; }, true );

    io::outputCMap( bdry, [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
        return sweep.mesh.points.at( indexingOrError( bdry, 0 )( v ) );
    }, "bdry.vtu" );

    const auto bdry_vert_ids = indexingOrError( bdry, 0 );
    const auto target_vert_ids = indexingOrError( target, 0 );

    std::array<topology::Vertex, 3> cut_vertices;
    iterateCellsWhile( bdry, 0, [&]( const topology::Vertex& v ) {
        if( bdry_vert_ids( v ) == cut_vert_ids.at( 0 ) ) cut_vertices.at( 0 ) = v;
        else if( bdry_vert_ids( v ) == cut_vert_ids.at( 1 ) ) cut_vertices.at( 1 ) = v;
        else if( bdry_vert_ids( v ) == cut_vert_ids.at( 2 ) ) cut_vertices.at( 2 ) = v;
        return true;
    } );

    const auto positions = [&]( const topology::Vertex& v ) {
        return sweep.mesh.points.at( bdry_vert_ids( v ) );
    };

    const auto cut1 = topology::shortestPath( target, positions, cut_vertices.at( 0 ), testEqualVertices( bdry_vert_ids, cut_vertices.at( 1 ) ) );
    topology::GlobalCellMarker cut_marker( bdry, 1 );
    for( const auto& e : cut1 ) cut_marker.mark( bdry, e );
    const auto cut2 = topology::shortestPath( target, [&]( const topology::Edge& e ){
        return cut_marker.isMarked( e ) ? std::numeric_limits<double>::max() : edgeLength( bdry, positions, e );
    }, cut_vertices.at( 1 ), testEqualVertices( bdry_vert_ids, cut_vertices.at( 2 ) ) );

    io::outputEdges( bdry, positions, util::concatenate( cut1, cut2 ), "level_set_cut.vtu" );

    std::set<topology::Cell> cuts;
    cuts.insert( cut1.begin(), cut1.end() );
    cuts.insert( cut2.begin(), cut2.end() );
    const topology::CutCombinatorialMap cut_cmap( target, cuts );
    const auto cutmap_vert_ids = indexingOrError( cut_cmap, 0 );

    const Eigen::MatrixX2d tutte = reparam::tutteOrbifoldEmbedding( cut_cmap, [&]( const topology::Vertex& v ) {
        return sweep.mesh.points.at( bdry_vert_ids( v ) );
    }, { cut_vertices.at( 0 ), cut_vertices.at( 1 ), cut_vertices.at( 2 ) }, reparam::Laplace2dEdgeWeights::InverseLength );

    io::outputCMap( cut_cmap, [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
        return Eigen::Vector3d( tutte( cutmap_vert_ids( v ), 0 ), tutte( cutmap_vert_ids( v ), 1 ), 0 );
    }, "tutte_orbi.vtu" );
}