#include <catch2/catch_test_macros.hpp>
#include <TriMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <SimplicialComplexTestCases.hpp>

using namespace topology;

TEST_CASE( "TriMeshSimplicialComplex on a 4-tri complex has correct number of cells" )
{
    const SimplicialComplex mesh = TriMeshTestCases::fourTriangleQuad();

    const TriMeshCombinatorialMap map( mesh );

    CHECK( cellCount( map, 2 ) == 4 );
    CHECK( cellCount( map, 1 ) == 8 );
    CHECK( cellCount( map, 0 ) == 5 );

    size_t n_boundary_edge = 0;
    iterateCellsWhile( map, 1, [&]( const Edge& e ) {
        if( boundaryAdjacent( map, e ) ) n_boundary_edge++;
        return true;
    } );
    CHECK( n_boundary_edge == 4 );

    size_t n_boundary_vert = 0;
    iterateCellsWhile( map, 0, [&]( const Vertex& v ) {
        if( boundaryAdjacent( map, v ) ) n_boundary_vert++;
        return true;
    } );
    CHECK( n_boundary_vert == 4 );
}

TEST_CASE( "TriMeshSimplicialComplex on a complex that is the boundary of a tet" )
{
    const SimplicialComplex mesh = TriMeshTestCases::tetBoundary();

    const TriMeshCombinatorialMap map( mesh );

    CHECK( cellCount( map, 2 ) == 4 );
    CHECK( cellCount( map, 1 ) == 6 );
    CHECK( cellCount( map, 0 ) == 4 );

    size_t n_boundary_edge = 0;
    iterateCellsWhile( map, 1, [&]( const Edge& e ) {
        if( boundaryAdjacent( map, e ) ) n_boundary_edge++;
        return true;
    } );
    CHECK( n_boundary_edge == 0 );

    size_t n_boundary_vert = 0;
    iterateCellsWhile( map, 0, [&]( const Vertex& v ) {
        if( boundaryAdjacent( map, v ) ) n_boundary_vert++;
        return true;
    } );
    CHECK( n_boundary_vert == 0 );
}