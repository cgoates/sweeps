#include <catch2/catch_test_macros.hpp>
#include <MultiPatchSplineFactory.hpp>
#include <TraceExtraction.hpp>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>
#include <algorithm>
#include <set>
#include <vector>

using namespace basis;
using namespace topology;

namespace
{
    using InternalConnectionsMap = MultiPatchCombinatorialMap::InternalConnectionsMap;
    using ConstituentSide = MultiPatchCombinatorialMap::ConstituentSide;
    using TPPermutation = MultiPatchCombinatorialMap::TPPermutation;

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

    InternalConnectionsMap twoPatchConnection( const ElementSide& first_side,
                                               const ElementSide& second_side,
                                               const TPPermutation permutation )
    {
        const ConstituentSide first{ 0, first_side.sideId() };
        const ConstituentSide second{ 1, second_side.sideId() };
        return { { first, { permutation, second } }, { second, { permutation, first } } };
    }

    std::shared_ptr<const TPSplineSpace> makePatch( const SmallVector<KnotVector, 3>& kvs,
                                                    const SmallVector<size_t, 3>& degrees )
    {
        return std::make_shared<const TPSplineSpace>( buildBSpline( kvs, degrees ) );
    }

    std::set<Eigen::Index> intersection( const std::set<Eigen::Index>& first,
                                         const std::set<Eigen::Index>& second )
    {
        std::set<Eigen::Index> out;
        std::set_intersection(
            first.begin(), first.end(), second.begin(), second.end(), std::inserter( out, out.end() ) );
        return out;
    }

    size_t intersectionSize( const std::set<Eigen::Index>& first, const std::set<Eigen::Index>& second )
    {
        return intersection( first, second ).size();
    }

    std::vector<Cell> topCellsOnSide( const TPCombinatorialMap& cmap, const ElementSide& side )
    {
        const SmallVector<std::shared_ptr<const CombinatorialMap1d>, 3> components =
            tensorProductComponentCMaps( cmap );

        util::IndexVec full_lengths;
        util::IndexVec side_lengths;
        for( size_t axis = 0; axis < cmap.dim(); axis++ )
        {
            const size_t length = cellCount( *components.at( axis ), components.at( axis )->dim() );
            full_lengths.push_back( length );
            if( axis != side.axis ) side_lengths.push_back( length );
        }

        std::vector<Cell> out;
        util::iterateTensorProduct( side_lengths, [&]( const util::IndexVec& side_indices ) {
            util::IndexVec element_indices;
            for( size_t axis = 0, side_axis = 0; axis < cmap.dim(); axis++ )
            {
                if( axis == side.axis )
                    element_indices.push_back( side.lower ? 0 : full_lengths.at( axis ) - 1 );
                else
                    element_indices.push_back( side_indices.at( side_axis++ ) );
            }
            out.push_back( topCellAt( cmap, element_indices ) );
        } );

        return out;
    }

    std::set<Eigen::Index> scalarGlobalSideIds( const MultiPatchSplineSpace& ss,
                                                const size_t patch_id,
                                                const ElementSide& side )
    {
        const TPSplineSpace& patch = *ss.subSpaces().at( patch_id );
        const auto& cmap = dynamic_cast<const TPCombinatorialMap&>( patch.basisComplex().parametricAtlas().cmap() );

        std::set<Eigen::Index> out;
        for( const Cell& cell : topCellsOnSide( cmap, side ) )
        {
            const TraceSideData trace = traceSideData( patch, cell, side );
            for( const FunctionId local_fid : trace.connectivity )
            {
                out.insert( ss.functionIdMap().at( patch_id ).at( local_fid.id() ).id() );
            }
        }
        return out;
    }

    size_t componentIndex( const VectorConformingTPSplineSpace& ss, const FunctionId local_fid )
    {
        Eigen::Index offset = 0;
        for( size_t component = 0; component < ss.scalarTPBases().size(); component++ )
        {
            const Eigen::Index next_offset =
                offset + static_cast<Eigen::Index>( ss.scalarTPBases().at( component )->numFunctions() );
            if( local_fid.id() < next_offset ) return component;
            offset = next_offset;
        }

        throw std::runtime_error( "Local vector function id is outside the component ranges." );
    }

    std::vector<std::set<Eigen::Index>> vectorGlobalSideIdsByComponent(
        const VectorConformingMultiPatchSplineSpace& ss,
        const size_t patch_id,
        const ElementSide& side )
    {
        const VectorConformingTPSplineSpace& patch = *ss.subSpaces().at( patch_id );
        const auto& cmap = dynamic_cast<const TPCombinatorialMap&>( patch.basisComplex().parametricAtlas().cmap() );

        std::vector<std::set<Eigen::Index>> out( patch.numVectorComponents() );
        for( const Cell& cell : topCellsOnSide( cmap, side ) )
        {
            const TraceSideData trace = traceSideData( patch, cell, side );
            for( const FunctionId local_fid : trace.connectivity )
            {
                const size_t component = componentIndex( patch, local_fid );
                out.at( component )
                    .insert( ss.functionIdMap().at( patch_id ).at( local_fid.id() ).first.id() );
            }
        }
        return out;
    }
}

TEST_CASE( "Phase 3 factories produce expected 2D multipatch DOF sharing" )
{
    const double ptol = 1e-10;
    const KnotVector kv( { 0, 0, 0, 1, 1, 1 }, ptol );
    const auto first_patch = makePatch( { kv, kv }, { 2, 2 } );
    const auto second_patch = makePatch( { kv, kv }, { 2, 2 } );

    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );
    const InternalConnectionsMap connections =
        twoPatchConnection( first_side, second_side, TPPermutation::Flip1d );

    const MultiPatchSplineSpace h1 =
        buildH1MultiPatchSplineSpace( { first_patch, second_patch }, connections );
    CHECK( h1.numFunctions() == 15 );
    CHECK( intersectionSize( scalarGlobalSideIds( h1, 0, first_side ),
                             scalarGlobalSideIds( h1, 1, second_side ) ) == 3 );

    const VectorConformingMultiPatchSplineSpace hdiv = buildHDivMultiPatchSplineSpace( h1 );
    const auto hdiv_first = vectorGlobalSideIdsByComponent( hdiv, 0, first_side );
    const auto hdiv_second = vectorGlobalSideIdsByComponent( hdiv, 1, second_side );
    REQUIRE( hdiv_first.size() == 2 );
    CHECK( hdiv.numFunctions() == 22 );
    CHECK( intersectionSize( hdiv_first.at( 0 ), hdiv_second.at( 0 ) ) == 2 );
    CHECK( intersectionSize( hdiv_first.at( 1 ), hdiv_second.at( 1 ) ) == 0 );

    const VectorConformingMultiPatchSplineSpace explicit_hdiv =
        buildVectorConformingMultiPatchSplineSpace( h1, ConformingType::Divergence );
    const auto explicit_hdiv_first = vectorGlobalSideIdsByComponent( explicit_hdiv, 0, first_side );
    const auto explicit_hdiv_second = vectorGlobalSideIdsByComponent( explicit_hdiv, 1, second_side );
    CHECK( explicit_hdiv.numFunctions() == hdiv.numFunctions() );
    CHECK( intersectionSize( explicit_hdiv_first.at( 0 ), explicit_hdiv_second.at( 0 ) ) == 2 );
    CHECK( intersectionSize( explicit_hdiv_first.at( 1 ), explicit_hdiv_second.at( 1 ) ) == 0 );

    const VectorConformingMultiPatchSplineSpace hcurl = buildHCurlMultiPatchSplineSpace( h1 );
    const auto hcurl_first = vectorGlobalSideIdsByComponent( hcurl, 0, first_side );
    const auto hcurl_second = vectorGlobalSideIdsByComponent( hcurl, 1, second_side );
    REQUIRE( hcurl_first.size() == 2 );
    CHECK( hcurl.numFunctions() == 22 );
    CHECK( intersectionSize( hcurl_first.at( 0 ), hcurl_second.at( 0 ) ) == 0 );
    CHECK( intersectionSize( hcurl_first.at( 1 ), hcurl_second.at( 1 ) ) == 2 );

    const VectorConformingMultiPatchSplineSpace explicit_hcurl =
        buildVectorConformingMultiPatchSplineSpace( h1, ConformingType::Curl );
    const auto explicit_hcurl_first = vectorGlobalSideIdsByComponent( explicit_hcurl, 0, first_side );
    const auto explicit_hcurl_second = vectorGlobalSideIdsByComponent( explicit_hcurl, 1, second_side );
    CHECK( explicit_hcurl.numFunctions() == hcurl.numFunctions() );
    CHECK( intersectionSize( explicit_hcurl_first.at( 0 ), explicit_hcurl_second.at( 0 ) ) == 0 );
    CHECK( intersectionSize( explicit_hcurl_first.at( 1 ), explicit_hcurl_second.at( 1 ) ) == 2 );

    const MultiPatchSplineSpace l2 = buildL2MultiPatchSplineSpace( hdiv );
    CHECK( l2.numFunctions() == 8 );
    CHECK( intersectionSize( scalarGlobalSideIds( l2, 0, first_side ),
                             scalarGlobalSideIds( l2, 1, second_side ) ) == 0 );
}

TEST_CASE( "Phase 3 factories handle multi-element 2D patch interfaces" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector kv_t_c0( { 0, 0, 0, 1, 1, 2, 2, 2 }, ptol );
    const auto first_patch = makePatch( { kv_s, kv_t_c0 }, { 2, 2 } );
    const auto second_patch = makePatch( { kv_s, kv_t_c0 }, { 2, 2 } );

    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );
    const InternalConnectionsMap connections =
        twoPatchConnection( first_side, second_side, TPPermutation::Flip1d );

    const MultiPatchSplineSpace h1 =
        buildH1MultiPatchSplineSpace( { first_patch, second_patch }, connections );
    CHECK( h1.numFunctions() == 25 );
    CHECK( intersectionSize( scalarGlobalSideIds( h1, 0, first_side ),
                             scalarGlobalSideIds( h1, 1, second_side ) ) == 5 );

    const auto& local_h1_cmap =
        dynamic_cast<const TPCombinatorialMap&>( h1.subSpaces().front()->basisComplex().parametricAtlas().cmap() );
    const TraceInterfaceElement local_h1_interface =
        traceInterfaceElement( *h1.subSpaces().front(), topCellAt( local_h1_cmap, { 0, 0 } ), ElementSide( 1, false ) );
    REQUIRE_FALSE( local_h1_interface.isBoundary() );
    REQUIRE( local_h1_interface.second.has_value() );

    std::set<Eigen::Index> local_h1_first_ids;
    for( const FunctionId fid : local_h1_interface.first.connectivity ) local_h1_first_ids.insert( fid.id() );
    std::set<Eigen::Index> local_h1_second_ids;
    for( const FunctionId fid : local_h1_interface.second->connectivity ) local_h1_second_ids.insert( fid.id() );
    CHECK( intersectionSize( local_h1_first_ids, local_h1_second_ids ) == 3 );

    const VectorConformingMultiPatchSplineSpace hdiv = buildHDivMultiPatchSplineSpace( h1 );
    const auto hdiv_first = vectorGlobalSideIdsByComponent( hdiv, 0, first_side );
    const auto hdiv_second = vectorGlobalSideIdsByComponent( hdiv, 1, second_side );
    REQUIRE( hdiv_first.size() == 2 );
    CHECK( hdiv.numFunctions() == 40 );
    CHECK( intersectionSize( hdiv_first.at( 0 ), hdiv_second.at( 0 ) ) == 4 );
    CHECK( intersectionSize( hdiv_first.at( 1 ), hdiv_second.at( 1 ) ) == 0 );

    const VectorConformingMultiPatchSplineSpace hcurl = buildHCurlMultiPatchSplineSpace( h1 );
    const auto hcurl_first = vectorGlobalSideIdsByComponent( hcurl, 0, first_side );
    const auto hcurl_second = vectorGlobalSideIdsByComponent( hcurl, 1, second_side );
    REQUIRE( hcurl_first.size() == 2 );
    CHECK( hcurl.numFunctions() == 40 );
    CHECK( intersectionSize( hcurl_first.at( 0 ), hcurl_second.at( 0 ) ) == 0 );
    CHECK( intersectionSize( hcurl_first.at( 1 ), hcurl_second.at( 1 ) ) == 4 );

    const MultiPatchSplineSpace l2 = buildL2MultiPatchSplineSpace( hdiv );
    CHECK( l2.numFunctions() == 16 );
    CHECK( intersectionSize( scalarGlobalSideIds( l2, 0, first_side ),
                             scalarGlobalSideIds( l2, 1, second_side ) ) == 0 );
}

TEST_CASE( "Phase 3 factories produce expected 3D multipatch face sharing" )
{
    const double ptol = 1e-10;
    const KnotVector kv( { 0, 0, 0, 1, 1, 1 }, ptol );
    const auto first_patch = makePatch( { kv, kv, kv }, { 2, 2, 2 } );
    const auto second_patch = makePatch( { kv, kv, kv }, { 2, 2, 2 } );

    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );
    const InternalConnectionsMap connections =
        twoPatchConnection( first_side, second_side, TPPermutation::ZeroToZero );

    const MultiPatchSplineSpace h1 =
        buildH1MultiPatchSplineSpace( { first_patch, second_patch }, connections );
    CHECK( h1.numFunctions() == 45 );
    CHECK( intersectionSize( scalarGlobalSideIds( h1, 0, first_side ),
                             scalarGlobalSideIds( h1, 1, second_side ) ) == 9 );

    const VectorConformingMultiPatchSplineSpace hdiv = buildHDivMultiPatchSplineSpace( h1 );
    const auto hdiv_first = vectorGlobalSideIdsByComponent( hdiv, 0, first_side );
    const auto hdiv_second = vectorGlobalSideIdsByComponent( hdiv, 1, second_side );
    REQUIRE( hdiv_first.size() == 3 );
    CHECK( hdiv.numFunctions() == 68 );
    CHECK( intersectionSize( hdiv_first.at( 0 ), hdiv_second.at( 0 ) ) == 4 );
    CHECK( intersectionSize( hdiv_first.at( 1 ), hdiv_second.at( 1 ) ) == 0 );
    CHECK( intersectionSize( hdiv_first.at( 2 ), hdiv_second.at( 2 ) ) == 0 );

    const VectorConformingMultiPatchSplineSpace hcurl = buildHCurlMultiPatchSplineSpace( h1 );
    const auto hcurl_first = vectorGlobalSideIdsByComponent( hcurl, 0, first_side );
    const auto hcurl_second = vectorGlobalSideIdsByComponent( hcurl, 1, second_side );
    REQUIRE( hcurl_first.size() == 3 );
    CHECK( hcurl.numFunctions() == 96 );
    CHECK( intersectionSize( hcurl_first.at( 0 ), hcurl_second.at( 0 ) ) == 0 );
    CHECK( intersectionSize( hcurl_first.at( 1 ), hcurl_second.at( 1 ) ) == 6 );
    CHECK( intersectionSize( hcurl_first.at( 2 ), hcurl_second.at( 2 ) ) == 6 );

    const MultiPatchSplineSpace l2 = buildL2MultiPatchSplineSpace( hdiv );
    CHECK( l2.numFunctions() == 16 );
    CHECK( intersectionSize( scalarGlobalSideIds( l2, 0, first_side ),
                             scalarGlobalSideIds( l2, 1, second_side ) ) == 0 );
}

TEST_CASE( "Phase 3 L2 factory keeps reduced C0 directions broken inside a patch" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s_c0( { 0, 0, 0, 1, 1, 2, 2, 2 }, ptol );
    const KnotVector kv_t( { 0, 0, 0, 1, 1, 1 }, ptol );
    const auto patch = makePatch( { kv_s_c0, kv_t }, { 2, 2 } );

    const MultiPatchSplineSpace h1 = buildH1MultiPatchSplineSpace( { patch }, InternalConnectionsMap{} );
    const VectorConformingMultiPatchSplineSpace hdiv = buildHDivMultiPatchSplineSpace( h1 );
    const MultiPatchSplineSpace l2 = buildL2MultiPatchSplineSpace( hdiv );

    REQUIRE( l2.subSpaces().size() == 1 );
    const TPSplineSpace& local_l2 = *l2.subSpaces().front();
    const auto& cmap = dynamic_cast<const TPCombinatorialMap&>( local_l2.basisComplex().parametricAtlas().cmap() );
    const TraceInterfaceElement iface =
        traceInterfaceElement( local_l2, topCellAt( cmap, { 0, 0 } ), ElementSide( 0, false ) );

    REQUIRE_FALSE( iface.isBoundary() );
    REQUIRE( iface.second.has_value() );

    std::set<Eigen::Index> first_ids;
    for( const FunctionId fid : iface.first.connectivity ) first_ids.insert( fid.id() );

    std::set<Eigen::Index> second_ids;
    for( const FunctionId fid : iface.second->connectivity ) second_ids.insert( fid.id() );

    CHECK( intersectionSize( first_ids, second_ids ) == 0 );
}
