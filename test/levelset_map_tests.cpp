#include <catch2/catch_test_macros.hpp>
#include <TriMeshCombinatorialMap.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <LevelSetCMap.hpp>
#include <SimplicialComplexTestCases.hpp>

using namespace topology;

TEST_CASE( "Various level set cmaps crossing 5 edges on a 4-tri complex" )
{
    const SimplicialComplex mesh = TriMeshTestCases::fourTriangleQuad();

    const TriMeshCombinatorialMap map( mesh );

    const Eigen::VectorXd func_values = ( Eigen::VectorXd( 5 ) << 0.4, 0, 0.9, 0.2, 1 ).finished();

    const auto vertex_ids = indexingOrError( map, 0 );

    const auto func = [&]( const Vertex& v ){
        return func_values( vertex_ids( v ) );
    };

    const auto check_phi1s_nonperiodic = []( const CombinatorialMap& level ) {
        bool found_forward_end = false;
        bool found_backward_end = false;
        iterateDartsWhile( level, [&]( const Dart& d ) {
            const auto phi1 = phi( level, 1, d );
            const auto phi_1 = phi( level, -1, d );
            if( found_backward_end ) CHECK( phi_1.has_value() );
            if( found_forward_end ) CHECK( phi1.has_value() );
            found_forward_end = found_forward_end or not phi1.has_value();
            found_backward_end = found_backward_end or not phi_1.has_value();
            return true;
        } );
    };

    {
        const LevelSetCMap level( map, func, 0.1 );
        CHECK( level.dim() == 1 );
        CHECK( cellCount( level, 1 ) == 2 );
        CHECK( cellCount( level, 0 ) == 2 );
        check_phi1s_nonperiodic( level );
    }
    {
        const LevelSetCMap level( map, func, 0.3 );
        CHECK( level.dim() == 1 );
        CHECK( cellCount( level, 1 ) == 3 );
        CHECK( cellCount( level, 0 ) == 3 );
        check_phi1s_nonperiodic( level );
    }
    {
        const LevelSetCMap level( map, func, 0.5 );
        CHECK( level.dim() == 1 );
        CHECK( cellCount( level, 1 ) == 4 );
        CHECK( cellCount( level, 0 ) == 4 );
        check_phi1s_nonperiodic( level );
    }
    {
        const LevelSetCMap level( map, func, 0.95 );
        CHECK( level.dim() == 1 );
        CHECK( cellCount( level, 1 ) == 4 );
        CHECK( cellCount( level, 0 ) == 4 );
        iterateDartsWhile( level, [&]( const Dart& d ) {
            CHECK( phi( level, 1, d ).has_value() );
            CHECK( phi( level, -1, d ).has_value() );
            return true;
        } );
    }
}

TEST_CASE( "Level sets on a two-tet complex" )
{
    const SimplicialComplex mesh = SweepInputTestCases::twoTets().mesh;
    const TetMeshCombinatorialMap map( mesh );

    const Eigen::VectorXd func_values = ( Eigen::VectorXd( 5 ) << 0, 1, 0, 0, 1 ).finished();
    const auto vertex_ids = indexingOrError( map, 0 );
    const auto func = [&]( const Vertex& v ) { return func_values( vertex_ids( v ) ); };

    {
        const LevelSetCMap level( map, func, 0.5 );
        CHECK( level.dim() == 2 );
        CHECK( cellCount( level, 2 ) == 2 );
        CHECK( cellCount( level, 1 ) == 6 );
        CHECK( cellCount( level, 0 ) == 5 );
        std::set<size_t> n_edges_in_faces;
        iterateCellsWhile( level, 2, [&]( const Face& f ) {
            size_t n_edges = 0;
            iterateAdjacentCells( level, f, 1, [&]( const auto& ){ n_edges++; return true; } );
            n_edges_in_faces.insert( n_edges );
            return true;
        } );
        CHECK( n_edges_in_faces == std::set<size_t>{ 3, 4 } );

        size_t n_darts_with_phi2 = 0;
        size_t n_darts = 0;
        iterateDartsWhile( level, [&]( const Dart& d ) {
            CHECK( phi( level, 1, d ).has_value() );
            CHECK( phi( level, -1, d ).has_value() );
            const auto noop = phi( level, {1, -1}, d );
            CHECK( noop.has_value() );
            if( noop.has_value() ) CHECK( noop.value() == d );
            if( phi( level, 2, d ).has_value() ) n_darts_with_phi2++;
            n_darts++;
            return true;
        } );
        CHECK( n_darts == 7 );
        CHECK( n_darts_with_phi2 == 2 );
    }
}