#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <CutCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <TriMeshCombinatorialMap.hpp>
#include <Logging.hpp>

using namespace topology;
TEST_CASE( "Cuts of a TP cmap" )
{
    const auto source = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto line = std::make_shared<const CombinatorialMap1d>( 4 );
    const TPCombinatorialMap tp( source, line );
    {
        const CutCombinatorialMap cut_cmap( tp, { Edge( 5 ), Edge( 17 ) } );
        CHECK( cellCount( tp, 0 ) + 2 == cellCount( cut_cmap, 0 ) );
        CHECK( cellCount( tp, 1 ) + 2 == cellCount( cut_cmap, 1 ) );
        CHECK( cellCount( tp, 2 ) == cellCount( cut_cmap, 2 ) );
    }
    {
        const CutCombinatorialMap cut_cmap( tp, { Edge( 17 ) } );
        CHECK( cellCount( tp, 0 ) == cellCount( cut_cmap, 0 ) );
        CHECK( cellCount( tp, 1 ) + 1 == cellCount( cut_cmap, 1 ) );
        CHECK( cellCount( tp, 2 ) == cellCount( cut_cmap, 2 ) );
    }
}

TEST_CASE( "Cuts of a cmap with vertex ids" )
{
    const SimplicialComplex sc = TriMeshTestCases::fourTriangleQuad();
    const TriMeshCombinatorialMap cmap( sc );

    std::set<Cell> cut_edges;
    iterateCellsWhile( cmap, 1, [&]( const Edge& e ) {
        if( not onBoundary( cmap, e.dart() ) )
        {
            cut_edges.insert( e );
            return false;
        }
        return true;
    } );

    const CutCombinatorialMap cut_cmap( cmap, cut_edges );
    CHECK( cellCount( cmap, 0 ) + 1 == cellCount( cut_cmap, 0 ) );
    CHECK( cellCount( cmap, 1 ) + 1 == cellCount( cut_cmap, 1 ) );
    CHECK( cellCount( cmap, 2 ) == cellCount( cut_cmap, 2 ) );

    std::set<size_t> actual_vert_ids;
    const auto cut_vert_ids = indexingOrError( cut_cmap, 0 );
    iterateCellsWhile( cut_cmap, 0, [&]( const Vertex& v ) {
        actual_vert_ids.insert( cut_vert_ids( v ) );
        return true;
    } );

    CHECK( actual_vert_ids == std::set<size_t>( { 0, 1, 2, 3, 4, 5 } ) );
}
