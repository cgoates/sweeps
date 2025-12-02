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

static constexpr bool LAPLACE_TESTS_OUTPUT_VTK = false;

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

    if( LAPLACE_TESTS_OUTPUT_VTK )
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

    if( LAPLACE_TESTS_OUTPUT_VTK )
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

    if( LAPLACE_TESTS_OUTPUT_VTK )
        io::outputEdgeChain( bdry, positions, util::concatenate( cut1, cut2 ), "level_set_cut.vtu" );

    std::set<topology::Cell> cuts;
    cuts.insert( cut1.begin(), cut1.end() );
    cuts.insert( cut2.begin(), cut2.end() );
    const topology::CutCombinatorialMap cut_cmap( target, cuts );
    const auto cutmap_vert_ids = indexingOrError( cut_cmap, 0 );

    const Eigen::MatrixX2d tutte = reparam::tutteOrbifoldEmbedding( cut_cmap, [&]( const topology::Vertex& v ) {
        return sweep.mesh.points.at( bdry_vert_ids( v ) );
    }, { cut_vertices.at( 0 ), cut_vertices.at( 1 ), cut_vertices.at( 2 ) }, reparam::Laplace2dEdgeWeights::InverseLength );

    if( LAPLACE_TESTS_OUTPUT_VTK )
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

    if( LAPLACE_TESTS_OUTPUT_VTK )
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

    if( LAPLACE_TESTS_OUTPUT_VTK )
        io::outputEdgeChain( bdry, positions, util::concatenate( cut1, cut2 ), "level_set_cut.vtu" );

    std::set<topology::Cell> cuts;
    cuts.insert( cut1.begin(), cut1.end() );
    cuts.insert( cut2.begin(), cut2.end() );
    const topology::CutCombinatorialMap cut_cmap( target, cuts );
    const auto cutmap_vert_ids = indexingOrError( cut_cmap, 0 );

    const Eigen::MatrixX2d tutte = reparam::tutteOrbifoldEmbedding( cut_cmap, [&]( const topology::Vertex& v ) {
        return sweep.mesh.points.at( bdry_vert_ids( v ) );
    }, { cut_vertices.at( 0 ), cut_vertices.at( 1 ), cut_vertices.at( 2 ) }, reparam::Laplace2dEdgeWeights::InverseLength );

    if( LAPLACE_TESTS_OUTPUT_VTK )
        io::outputCMap( cut_cmap, [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
            return Eigen::Vector3d( tutte( cutmap_vert_ids( v ), 0 ), tutte( cutmap_vert_ids( v ), 1 ), 0 );
        }, "tutte_orbi.vtu" );
}

SimplicialComplex annulusExample()
{
    std::vector<Eigen::Vector3d> V;
    V.reserve( 15 );

    auto add_ring = [&]( double r ) {
        for( int i = 0; i < 5; ++i )
        {
            double theta = 2.0 * M_PI * i / 5.0; // evenly spaced pentagon
            V.emplace_back( r * std::cos( theta ), r * std::sin( theta ), 0.0 );
        }
    };

    add_ring( 1.0 ); // outer ring: indices 0–4
    add_ring( 0.7 ); // middle ring: indices 5–9
    add_ring( 0.4 ); // inner ring: indices 10–14

    V.at( 5 ) = Eigen::Vector3d( 0.7, 0.2, 0.0 );   // perturb middle ring a bit

    std::vector<Simplex> F;
    F.reserve( 20 );

    auto add_strip = [&]( size_t A0, size_t B0 ) {
        // A = outer ring start index
        // B = inner ring start index
        for( int i = 0; i < 5; ++i )
        {
            size_t ai = A0 + i;
            size_t an = A0 + ( i + 1 ) % 5;
            size_t bi = B0 + i;
            size_t bn = B0 + ( i + 1 ) % 5;

            F.push_back( { ai, an, bi } );
            F.push_back( { an, bn, bi } );
        }
    };

    // Outer ↔ middle
    add_strip( 0, 5 );

    // Middle ↔ inner
    add_strip( 5, 10 );

    return SimplicialComplex{ F, V };
}

bool isPositivelyOriented( const Triangle<2>& tri )
{
    const Eigen::Vector2d v0 = tri.v2 - tri.v1;
    const Eigen::Vector2d v1 = tri.v3 - tri.v1;
    return ( v0( 0 ) * v1( 1 ) - v0( 1 ) * v1( 0 ) ) > 0;
}

TEST_CASE( "Cut tutte embedding" )
{
    const SimplicialComplex mesh = annulusExample();

    const topology::TriMeshCombinatorialMap map( mesh );
    const topology::CutCombinatorialMap cut_map( map, { topology::Edge( 27 ), topology::Edge( 30 ), topology::Edge( 31 ) } );

    const auto sol = reparam::cutTutteEmbedding( cut_map,
                                                   [&]( const topology::Vertex& v ) {
                                                       const auto vids = indexingOrError( map, 0 );
                                                       return mesh.points.at( vids( v ) );
                                                   },
                                                   { { topology::Vertex( 27 ),topology::Vertex( 58 ) } },
                                                   reparam::Laplace2dEdgeWeights::InverseLength ).first;

    const auto cut_vert_ids = indexingOrError( cut_map, 0 );
    iterateCellsWhile( map, 2, [&]( const topology::Face& f ) {
        const auto tri = triangleOfFace<2>( cut_map, [&]( const auto& v ) -> Eigen::Vector2d { return sol.row( cut_vert_ids( v ) ); }, f );
        CHECK( not isPositivelyOriented( tri ) ); // Doesn't matter if it's positive, just that they're consistent
        return true;
    } );

    if( LAPLACE_TESTS_OUTPUT_VTK )
    {
        SimplicialComplex cut_mesh;
        cut_mesh.simplices = std::vector<Simplex>( mesh.simplices.size(), Simplex(0) );
        cut_mesh.points = std::vector<Eigen::Vector3d>( cellCount( cut_map, 0 ) );
        const auto vert_ids = indexingOrError( map, 0 );
        const auto cut_face_ids = indexingOrError( cut_map, 2 );
        iterateCellsWhile( cut_map, 0, [&]( const topology::Vertex& v ) {
            cut_mesh.points.at( cut_vert_ids( v ) ) = mesh.points.at( vert_ids( v ) );
            return true;
        } );
        iterateCellsWhile( map, 2, [&]( const topology::Face& f ) {
            cut_mesh.simplices.at( cut_face_ids( f ) ) = Simplex(
                cut_vert_ids( topology::Vertex( f.dart() ) ),
                cut_vert_ids( topology::Vertex( phi( map, 1, f.dart() ).value() ) ),
                cut_vert_ids( topology::Vertex( phi( map, -1, f.dart() ).value() ) ) );
            return true;
        } );

        Eigen::MatrixX3d sol_3d( sol.rows(), 3 );
        sol_3d.col( 0 ) = sol.col( 0 );
        sol_3d.col( 1 ) = sol.col( 1 );
        sol_3d.col( 2 ).setZero();

        io::VTKOutputObject out( cut_mesh );
        out.addVertexField( "tutte_embedding", sol_3d );
        io::outputSimplicialFieldToVTK( out, "input_mesh.vtu" );
    }
}

TEST_CASE( "Cut tutte embedding multiple cuts" )
{
    const SimplicialComplex mesh = []() {
        // Load the mesh from the OBJ file
        const auto mesh = io::loadOBJFile( SRC_HOME "/test/data/doubleAnnulus.obj" );
        SimplicialComplex sc;
        for( const auto& conn : mesh.first )
        {
            sc.simplices.push_back( Simplex( conn.at( 0 ), conn.at( 1 ), conn.at( 2 ) ) );
        }
        sc.points = mesh.second;
        return sc;
    }();

    // First cut: vertex indices 55, 228, 220, 12
    // Second cut: vertex indices 14, 113, 210, 108, 8

    std::vector<std::pair<topology::Vertex, topology::Vertex>> cut_extremities{
        { topology::Vertex( 999990 ),topology::Vertex( 999990 ) },
        { topology::Vertex( 999990 ),topology::Vertex( 999990 ) }
    };

    const topology::TriMeshCombinatorialMap map( mesh );
    const std::set<topology::Cell> cuts = [&](){
        const std::set<int> first_vids = { 55, 228, 220, 12 };
        const std::set<int> second_vids = { 14, 113, 210, 108, 8 };
        std::set<topology::Cell> result;
        const auto vids = indexingOrError( map, 0 );
        iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
            const topology::Vertex v1 = topology::Vertex( e.dart() );
            const topology::Vertex v2 = topology::Vertex( phi( map, 1, e.dart() ).value() );
            if( first_vids.count( vids( v1 ) ) and first_vids.count( vids( v2 ) ) )
            {
                result.insert( e );
            }
            if( second_vids.count( vids( v1 ) ) and second_vids.count( vids( v2 ) ) )
            {
                result.insert( e );
            }
            if( ( vids( v1 ) == 55 and vids( v2 ) == 228) or( vids( v1 ) == 228 and vids( v2 ) == 55 ) )
            {
                cut_extremities.at( 0 ).first = v1;
            }
            if( ( vids( v1 ) == 220 and vids( v2 ) == 12) or( vids( v1 ) == 12 and vids( v2 ) == 220 ) )
            {
                cut_extremities.at( 0 ).second = v2;
            }
            if( ( vids( v1 ) == 14 and vids( v2 ) == 113) or( vids( v1 ) == 113 and vids( v2 ) == 14 ) )
            {
                cut_extremities.at( 1 ).first = v1;
            }
            if( vids( v1 ) == 8 and vids( v2 ) == 108 )
            {
                // Find the dart furthest along the phi 2,1 cycle from v1
                topology::Vertex current_v = v1;
                do
                {
                    const auto next_dart_opt = phi( map, {2,1}, current_v.dart() );
                    if( not next_dart_opt.has_value() ) break;
                    current_v = topology::Vertex( next_dart_opt.value() );
                } while( vids( current_v ) != vids( v1 ) );
                cut_extremities.at( 1 ).second = current_v;
            }
            return true;
        } );
        return result;
    }();
    const topology::CutCombinatorialMap cut_map( map, cuts );

    const auto sol = reparam::cutTutteEmbedding( cut_map,
                                                   [&]( const topology::Vertex& v ) {
                                                       const auto vids = indexingOrError( map, 0 );
                                                       return mesh.points.at( vids( v ) );
                                                   },
                                                   cut_extremities,
                                                   reparam::Laplace2dEdgeWeights::InverseLength ).first;

    const auto cut_vert_ids = indexingOrError( cut_map, 0 );
    iterateCellsWhile( map, 2, [&]( const topology::Face& f ) {
        const auto tri = triangleOfFace<2>(
            cut_map, [&]( const auto& v ) -> Eigen::Vector2d { return sol.row( cut_vert_ids( v ) ); }, f );
        CHECK( not isPositivelyOriented( tri ) ); // Doesn't matter if it's positive, just that they're consistent
        return true;
    } );

    if( LAPLACE_TESTS_OUTPUT_VTK )
    {
        SimplicialComplex cut_mesh;
        cut_mesh.simplices = std::vector<Simplex>( mesh.simplices.size(), Simplex(0) );
        cut_mesh.points = std::vector<Eigen::Vector3d>( cellCount( cut_map, 0 ) );
        const auto vert_ids = indexingOrError( map, 0 );
        const auto cut_vert_ids = indexingOrError( cut_map, 0 );
        const auto cut_face_ids = indexingOrError( cut_map, 2 );
        iterateCellsWhile( cut_map, 0, [&]( const topology::Vertex& v ) {
            cut_mesh.points.at( cut_vert_ids( v ) ) = mesh.points.at( vert_ids( v ) );
            return true;
        } );
        iterateCellsWhile( map, 2, [&]( const topology::Face& f ) {
            cut_mesh.simplices.at( cut_face_ids( f ) ) = Simplex(
                cut_vert_ids( topology::Vertex( f.dart() ) ),
                cut_vert_ids( topology::Vertex( phi( map, 1, f.dart() ).value() ) ),
                cut_vert_ids( topology::Vertex( phi( map, -1, f.dart() ).value() ) ) );
            return true;
        } );

        Eigen::MatrixX3d sol_3d( sol.rows(), 3 );
        sol_3d.col( 0 ) = sol.col( 0 );
        sol_3d.col( 1 ) = sol.col( 1 );
        sol_3d.col( 2 ).setZero();

        io::VTKOutputObject out( cut_mesh );
        out.addVertexField( "tutte_embedding", sol_3d );
        io::outputSimplicialFieldToVTK( out, "input_mesh.vtu" );
    }
}
