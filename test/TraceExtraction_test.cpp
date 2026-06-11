#include <catch2/catch_test_macros.hpp>
#include <TraceExtraction.hpp>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>
#include <TPSplineSpace.hpp>
#include <VectorConformingBasisComplex.hpp>
#include <VectorConformingTPSplineSpace.hpp>
#include <algorithm>
#include <set>

using namespace basis;
using namespace topology;

namespace
{
    Cell topCellAt( const TPCombinatorialMap& cmap, const util::IndexVec& element_indices )
    {
        FullyUnflattenedDart unflat;
        for( const size_t element_index : element_indices )
        {
            unflat.unflat_darts.push_back( Dart( element_index ) );
        }
        unflat.dart_pos =
            SmallVector<TPCombinatorialMap::TPDartPos, 2>(
                element_indices.size() - 1, TPCombinatorialMap::TPDartPos::DartPos0 );
        return Cell( flattenFull( cmap, unflat ), cmap.dim() );
    }

    std::set<Eigen::Index> idSet( const std::vector<FunctionId>& conn )
    {
        std::set<Eigen::Index> out;
        std::transform( conn.begin(), conn.end(), std::inserter( out, out.end() ), []( const FunctionId& fid ) {
            return fid.id();
        } );
        return out;
    }

    size_t intersectionSize( const std::vector<FunctionId>& a, const std::vector<FunctionId>& b )
    {
        const std::set<Eigen::Index> a_ids = idSet( a );
        const std::set<Eigen::Index> b_ids = idSet( b );
        std::vector<Eigen::Index> intersection;
        std::set_intersection(
            a_ids.begin(), a_ids.end(), b_ids.begin(), b_ids.end(), std::back_inserter( intersection ) );
        return intersection.size();
    }

    void checkPartitionOfUnity( const Eigen::MatrixXd& extraction )
    {
        CHECK( ( extraction.colwise().sum().transpose() - Eigen::VectorXd::Ones( extraction.cols() ) ).norm() <
               1e-14 );
    }

    void checkSharedTraceRowsAgree( const TraceSideData& first, const TraceSideData& second, const double tol )
    {
        size_t num_shared = 0;
        for( size_t first_row = 0; first_row < first.connectivity.size(); first_row++ )
        {
            const auto second_it =
                std::find( second.connectivity.begin(), second.connectivity.end(), first.connectivity.at( first_row ) );
            if( second_it == second.connectivity.end() ) continue;

            const size_t second_row = std::distance( second.connectivity.begin(), second_it );
            CHECK( ( first.extraction.row( first_row ) - second.extraction.row( second_row ) ).norm() < tol );
            num_shared++;
        }

        CHECK( num_shared == intersectionSize( first.connectivity, second.connectivity ) );
    }
}

TEST_CASE( "ElementSide follows multipatch side id convention" )
{
    CHECK( ElementSide( 0, false ).sideId() == 0 );
    CHECK( ElementSide( 0, true ).sideId() == 1 );
    CHECK( ElementSide( 1, false ).sideId() == 2 );
    CHECK( ElementSide( 1, true ).sideId() == 3 );

    CHECK( elementSideFromId( 4 ) == ElementSide( 2, false ) );
    CHECK( elementSideFromId( 5 ) == ElementSide( 2, true ) );
    CHECK( ElementSide( 2, true ).opposite() == ElementSide( 2, false ) );
}

TEST_CASE( "Trace column indices select scalar tensor-product Bernstein side columns" )
{
    const ParentBasis surface_pb{ param::cubeDomain( 2 ), { bernsteinBasis( 2 ), bernsteinBasis( 1 ) } };

    CHECK( traceDegrees( degrees( surface_pb ), ElementSide( 0, true ) ) == SmallVector<size_t, 3>{ 1 } );
    CHECK( traceColumnIndices( surface_pb, ElementSide( 0, true ) ) == std::vector<Eigen::Index>{ 0, 3 } );
    CHECK( traceColumnIndices( surface_pb, ElementSide( 0, false ) ) == std::vector<Eigen::Index>{ 2, 5 } );
    CHECK( traceColumnIndices( surface_pb, ElementSide( 1, true ) ) == std::vector<Eigen::Index>{ 0, 1, 2 } );
    CHECK( traceColumnIndices( surface_pb, ElementSide( 1, false ) ) == std::vector<Eigen::Index>{ 3, 4, 5 } );

    const ParentBasis volume_pb{
        param::cubeDomain( 3 ), { bernsteinBasis( 1 ), bernsteinBasis( 2 ), bernsteinBasis( 1 ) } };

    CHECK( traceDegrees( degrees( volume_pb ), ElementSide( 0, true ) ) ==
           SmallVector<size_t, 3>{ 2, 1 } );
    CHECK( traceDegrees( degrees( volume_pb ), ElementSide( 1, true ) ) ==
           SmallVector<size_t, 3>{ 1, 1 } );
    CHECK( traceDegrees( degrees( volume_pb ), ElementSide( 2, false ) ) ==
           SmallVector<size_t, 3>{ 1, 2 } );
    CHECK( traceColumnIndices( volume_pb, ElementSide( 1, true ) ) ==
           std::vector<Eigen::Index>{ 0, 1, 6, 7 } );
    CHECK( traceColumnIndices( volume_pb, ElementSide( 1, false ) ) ==
           std::vector<Eigen::Index>{ 4, 5, 10, 11 } );
    CHECK( traceColumnIndices( volume_pb, ElementSide( 2, false ) ) ==
           std::vector<Eigen::Index>{ 6, 7, 8, 9, 10, 11 } );
}

TEST_CASE( "Trace column indices include vector-valued component blocks" )
{
    const ParentBasis div_pb = divConformingBernsteinCube( 2, 2 );
    CHECK( traceColumnIndices( div_pb, ElementSide( 0, true ) ) ==
           std::vector<Eigen::Index>{ 0, 3, 6, 8, 10 } );
    CHECK( traceColumnIndices( div_pb, ElementSide( 1, false ) ) ==
           std::vector<Eigen::Index>{ 3, 4, 5, 10, 11 } );

    const ParentBasis curl_pb = curlConformingBernsteinCube( 2, 2 );
    CHECK( traceColumnIndices( curl_pb, ElementSide( 0, true ) ) ==
           std::vector<Eigen::Index>{ 0, 2, 4, 6, 9 } );
    CHECK( traceColumnIndices( curl_pb, ElementSide( 1, false ) ) ==
           std::vector<Eigen::Index>{ 4, 5, 9, 10, 11 } );
}

TEST_CASE( "Trace side data prunes rows to one-sided side-supported functions" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s_c0( { 0, 0, 0, 1, 1, 2, 2, 2 }, ptol );
    const KnotVector kv_s_broken( { 0, 0, 0, 1, 1, 1, 2, 2, 2 }, ptol );
    const KnotVector kv_t( { 0, 0, 1, 1 }, ptol );

    SECTION( "C0 scalar traces on an interior side share the same trace DOFs" )
    {
        const TPSplineSpace ss = buildBSpline( { kv_s_c0, kv_t }, { 2, 1 } );
        const auto& cmap = dynamic_cast<const TPCombinatorialMap&>( ss.basisComplex().parametricAtlas().cmap() );
        const Cell first_elem = topCellAt( cmap, { 0, 0 } );

        const TraceInterfaceElement iface = traceInterfaceElement( ss, first_elem, ElementSide( 0, false ) );
        REQUIRE_FALSE( iface.isBoundary() );
        REQUIRE( iface.second.has_value() );
        CHECK( iface.first.connectivity.size() == 2 );
        CHECK( iface.second->connectivity.size() == 2 );
        CHECK( idSet( iface.first.connectivity ) == idSet( iface.second->connectivity ) );
        CHECK( iface.first.extraction.rows() == 2 );
        CHECK( iface.first.extraction.cols() == 2 );
        CHECK( iface.second->extraction.rows() == 2 );
        CHECK( iface.second->extraction.cols() == 2 );
        checkPartitionOfUnity( iface.first.extraction );
        checkPartitionOfUnity( iface.second->extraction );
        checkSharedTraceRowsAgree( iface.first, iface.second.value(), 1e-14 );

        CHECK( traceInterfaceElement( ss, first_elem, ElementSide( 0, true ) ).isBoundary() );
    }

    SECTION( "Broken scalar traces on an interior side keep independent one-sided DOFs" )
    {
        const TPSplineSpace ss = buildBSpline( { kv_s_broken, kv_t }, { 2, 1 } );
        const auto& cmap = dynamic_cast<const TPCombinatorialMap&>( ss.basisComplex().parametricAtlas().cmap() );
        const Cell first_elem = topCellAt( cmap, { 0, 0 } );

        const TraceInterfaceElement iface = traceInterfaceElement( ss, first_elem, ElementSide( 0, false ) );
        REQUIRE_FALSE( iface.isBoundary() );
        REQUIRE( iface.second.has_value() );
        CHECK( iface.first.connectivity.size() == 2 );
        CHECK( iface.second->connectivity.size() == 2 );
        CHECK( intersectionSize( iface.first.connectivity, iface.second->connectivity ) == 0 );
        checkPartitionOfUnity( iface.first.extraction );
        checkPartitionOfUnity( iface.second->extraction );
    }
}

TEST_CASE( "Trace side data supports vector-conforming tensor-product splines" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector kv_t( { 0, 0, 0, 1, 1, 1 }, ptol );
    const TPSplineSpace h1 = buildBSpline( { kv_s, kv_t }, { 2, 2 } );
    const auto hdiv_bc =
        std::make_shared<const VectorConformingBasisComplex>( h1.basisComplexPtr(), ConformingType::Divergence );
    const VectorConformingTPSplineSpace hdiv( hdiv_bc, h1 );

    const auto& cmap = dynamic_cast<const TPCombinatorialMap&>( h1.basisComplex().parametricAtlas().cmap() );
    const TraceSideData side = traceSideData( hdiv, topCellAt( cmap, { 0, 0 } ), ElementSide( 0, true ) );

    CHECK( side.connectivity.size() == 5 );
    CHECK( side.element_row_indices.size() == 5 );
    CHECK( side.extraction.rows() == 5 );
    CHECK( side.extraction.cols() == 5 );
    checkPartitionOfUnity( side.extraction );
}

TEST_CASE( "Explicit trace interface records support patch-local element sides" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector kv_t( { 0, 0, 1, 1 }, ptol );
    const TPSplineSpace first_patch = buildBSpline( { kv_s, kv_t }, { 2, 1 } );
    const TPSplineSpace second_patch = buildBSpline( { kv_s, kv_t }, { 2, 1 } );

    const auto& first_cmap =
        dynamic_cast<const TPCombinatorialMap&>( first_patch.basisComplex().parametricAtlas().cmap() );
    const auto& second_cmap =
        dynamic_cast<const TPCombinatorialMap&>( second_patch.basisComplex().parametricAtlas().cmap() );

    const TraceInterfaceElement iface = traceInterfaceElement( first_patch,
                                                               topCellAt( first_cmap, { 0, 0 } ),
                                                               ElementSide( 0, false ),
                                                               second_patch,
                                                               topCellAt( second_cmap, { 0, 0 } ),
                                                               ElementSide( 0, true ) );

    REQUIRE_FALSE( iface.isBoundary() );
    REQUIRE( iface.second.has_value() );
    CHECK( iface.first.connectivity.size() == 2 );
    CHECK( iface.second->connectivity.size() == 2 );
    CHECK( iface.first.extraction.rows() == 2 );
    CHECK( iface.first.extraction.cols() == 2 );
    CHECK( iface.second->extraction.rows() == 2 );
    CHECK( iface.second->extraction.cols() == 2 );
    checkPartitionOfUnity( iface.first.extraction );
    checkPartitionOfUnity( iface.second->extraction );
}

TEST_CASE( "Adjacent element side lookup distinguishes interior and boundary sides in 3D" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s( { 0, 0, 1, 2, 2 }, ptol );
    const KnotVector kv_t( { 0, 0, 1, 1 }, ptol );
    const KnotVector kv_u( { 0, 0, 1, 1 }, ptol );
    const TPSplineSpace ss = buildBSpline( { kv_s, kv_t, kv_u }, { 1, 1, 1 } );

    const auto& cmap = dynamic_cast<const TPCombinatorialMap&>( ss.basisComplex().parametricAtlas().cmap() );
    const Cell first_elem = topCellAt( cmap, { 0, 0, 0 } );
    const Cell second_elem = topCellAt( cmap, { 1, 0, 0 } );

    const auto neighbor = adjacentElementSide( cmap, first_elem, ElementSide( 0, false ) );
    REQUIRE( neighbor.has_value() );
    CHECK( neighbor->first.dart() == second_elem.dart() );
    CHECK( neighbor->second == ElementSide( 0, true ) );

    CHECK_FALSE( adjacentElementSide( cmap, first_elem, ElementSide( 0, true ) ).has_value() );
    CHECK_FALSE( adjacentElementSide( cmap, second_elem, ElementSide( 0, false ) ).has_value() );
}
