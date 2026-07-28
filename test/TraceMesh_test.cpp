#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <HierarchicalMultiPatchSplineSpace.hpp>
#include <MultiPatchSplineFactory.hpp>
#include <TraceMesh.hpp>
#include <VectorConformingHierarchicalMultiPatchSplineSpace.hpp>
#include <algorithm>
#include <set>

using namespace basis;
using namespace topology;

namespace
{
    using InternalConnectionsMap = MultiPatchCombinatorialMap::InternalConnectionsMap;
    using ConstituentSide = MultiPatchCombinatorialMap::ConstituentSide;
    using TPPermutation = MultiPatchCombinatorialMap::TPPermutation;

    Cell topCellAt( const TPCombinatorialMap& cmap, const std::vector<size_t>& element_indices )
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

    Cell patchTopCellAt( const MultiPatchCombinatorialMap& cmap,
                         const size_t patch_id,
                         const std::vector<size_t>& element_indices )
    {
        const Cell local = topCellAt( *cmap.constituents().at( patch_id ), element_indices );
        return Cell( cmap.toGlobalDart( patch_id, local.dart() ), cmap.dim() );
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

    MultiPatchSplineSpace makeTwoPatchSurface()
    {
        const double ptol = 1e-10;
        const KnotVector kv_s( { 0, 0, 0, 1, 1, 1 }, ptol );
        const KnotVector kv_t_c0( { 0, 0, 0, 1, 1, 2, 2, 2 }, ptol );
        const auto patch = makePatch( { kv_s, kv_t_c0 }, { 2, 2 } );
        return buildH1MultiPatchSplineSpace(
            { patch, patch },
            twoPatchConnection( ElementSide( 0, false ), ElementSide( 0, true ), TPPermutation::Flip1d ) );
    }

    MultiPatchSplineSpace makeTwoPatchVolume()
    {
        const double ptol = 1e-10;
        const KnotVector kv( { 0, 0, 0, 1, 1, 1 }, ptol );
        const auto patch = makePatch( { kv, kv, kv }, { 2, 2, 2 } );
        return buildH1MultiPatchSplineSpace(
            { patch, patch },
            twoPatchConnection( ElementSide( 0, false ), ElementSide( 0, true ), TPPermutation::ZeroToZero ) );
    }

    HierarchicalMultiPatchSplineSpace makeMixedLevelTwoPatchHierarchy()
    {
        const double ptol = 1e-10;
        const KnotVector coarse_s( { 0, 0, 0, 1, 2, 2, 2 }, ptol );
        const KnotVector fine_s( { 0, 0, 0, 0.5, 1, 1.5, 2, 2, 2 }, ptol );
        const KnotVector t( { 0, 0, 0, 1, 1, 1 }, ptol );

        const auto coarse_patch = makePatch( { coarse_s, t }, { 2, 2 } );
        const auto fine_patch = makePatch( { fine_s, t }, { 2, 2 } );
        const InternalConnectionsMap connections =
            twoPatchConnection( ElementSide( 0, false ), ElementSide( 0, true ), TPPermutation::Flip1d );

        auto coarse_h1 = std::make_shared<const MultiPatchSplineSpace>(
            buildH1MultiPatchSplineSpace( { coarse_patch, coarse_patch }, connections ) );
        auto fine_h1 = std::make_shared<const MultiPatchSplineSpace>(
            buildH1MultiPatchSplineSpace( { fine_patch, fine_patch }, connections ) );

        const auto& coarse_cmap =
            dynamic_cast<const MultiPatchCombinatorialMap&>(
                coarse_h1->basisComplex().parametricAtlas().cmap() );
        const auto& fine_cmap =
            dynamic_cast<const MultiPatchCombinatorialMap&>(
                fine_h1->basisComplex().parametricAtlas().cmap() );

        return buildHierarchicalSplineSpace(
            { coarse_h1, fine_h1 },
            {
                {
                    patchTopCellAt( coarse_cmap, 0, { 0, 0 } ),
                    patchTopCellAt( coarse_cmap, 1, { 1, 0 } )
                },
                {
                    patchTopCellAt( fine_cmap, 0, { 2, 0 } ),
                    patchTopCellAt( fine_cmap, 0, { 3, 0 } ),
                    patchTopCellAt( fine_cmap, 1, { 0, 0 } ),
                    patchTopCellAt( fine_cmap, 1, { 1, 0 } )
                }
            } );
    }

    size_t sharedTraceIdCount( const TraceInterfaceElement& iface )
    {
        REQUIRE( iface.second.has_value() );
        std::set<FunctionId> first( iface.first.connectivity.begin(), iface.first.connectivity.end() );
        std::set<FunctionId> second( iface.second->connectivity.begin(), iface.second->connectivity.end() );
        std::vector<FunctionId> shared;
        std::set_intersection( first.begin(), first.end(), second.begin(), second.end(), std::back_inserter( shared ) );
        return shared.size();
    }

    size_t aggregateSharedTraceIdCount( const SplineSpace& ss,
                                        const std::vector<TraceMeshInterface>& interfaces )
    {
        std::set<FunctionId> first_ids;
        std::set<FunctionId> second_ids;
        for( const TraceMeshInterface& iface : interfaces )
        {
            const TraceInterfaceElement trace = traceInterfaceElement( ss, iface );
            REQUIRE( trace.second.has_value() );
            first_ids.insert( trace.first.connectivity.begin(), trace.first.connectivity.end() );
            second_ids.insert( trace.second->connectivity.begin(), trace.second->connectivity.end() );
        }

        std::vector<FunctionId> shared;
        std::set_intersection( first_ids.begin(),
                               first_ids.end(),
                               second_ids.begin(),
                               second_ids.end(),
                               std::back_inserter( shared ) );
        return shared.size();
    }
}

TEST_CASE( "Phase 5A trace side coordinate permutations match multipatch conventions" )
{
    CHECK( permuteTraceSideIndex( { 0 }, { 3 }, TPPermutation::Flip1d ) == std::vector<size_t>{ 2 } );
    CHECK_THAT( permuteTraceSidePoint( Eigen::Vector<double, 1>( 0.25 ), TPPermutation::Flip1d )( 0 ),
                Catch::Matchers::WithinAbs( 0.75, 1e-14 ) );

    CHECK( permuteTraceSideIndex( { 0, 3 }, { 4, 5 }, TPPermutation::ZeroToZero ) ==
           std::vector<size_t>{ 3, 3 } );
    CHECK( permuteTraceSideIndex( { 0, 3 }, { 4, 5 }, TPPermutation::ZeroToOne ) ==
           std::vector<size_t>{ 1, 3 } );
    CHECK( permuteTraceSideIndex( { 0, 3 }, { 4, 5 }, TPPermutation::ZeroToTwo ) ==
           std::vector<size_t>{ 0, 1 } );
    CHECK( permuteTraceSideIndex( { 0, 3 }, { 4, 5 }, TPPermutation::ZeroToThree ) ==
           std::vector<size_t>{ 3, 0 } );

    const Eigen::Vector2d point( 0.2, 0.7 );
    CHECK( ( permuteTraceSidePoint( point, TPPermutation::ZeroToZero ) - Eigen::Vector2d( 0.8, 0.7 ) ).norm() <
           1e-14 );
    CHECK( ( permuteTraceSidePoint( point, TPPermutation::ZeroToOne ) - Eigen::Vector2d( 0.3, 0.8 ) ).norm() <
           1e-14 );
    CHECK( ( permuteTraceSidePoint( point, TPPermutation::ZeroToTwo ) - Eigen::Vector2d( 0.2, 0.3 ) ).norm() <
           1e-14 );
    CHECK( ( permuteTraceSidePoint( point, TPPermutation::ZeroToThree ) - Eigen::Vector2d( 0.7, 0.2 ) ).norm() <
           1e-14 );
}

TEST_CASE( "Phase 5A trace mesh enumerates simple 3D patch interfaces" )
{
    const MultiPatchSplineSpace h1 = makeTwoPatchVolume();

    const std::vector<TraceMeshInterface> boundaries = boundaryTraceMeshInterfaces( h1 );
    const std::vector<TraceMeshInterface> patch_interfaces = patchTraceMeshInterfaces( h1 );
    const std::vector<TraceMeshInterface> interiors = interiorTraceMeshInterfaces( h1 );
    const std::vector<TraceMeshInterface> all =
        traceMeshInterfaces( h1, true, true, true );

    CHECK( boundaries.size() == 10 );
    REQUIRE( patch_interfaces.size() == 1 );
    CHECK( interiors.empty() );
    CHECK( all.size() == boundaries.size() + patch_interfaces.size() + interiors.size() );

    const TraceMeshInterface& iface = patch_interfaces.front();
    CHECK( iface.type == TraceMeshEntityType::PatchInterface );
    REQUIRE_FALSE( iface.isBoundary() );
    CHECK( iface.first.patch_id == 0 );
    CHECK( iface.second->patch_id == 1 );
    CHECK( iface.first.side == ElementSide( 0, false ) );
    CHECK( iface.second->side == ElementSide( 0, true ) );
    REQUIRE( iface.permutation.has_value() );
    CHECK( iface.permutation.value() == TPPermutation::ZeroToZero );
    CHECK( sharedTraceIdCount( traceInterfaceElement( h1, iface ) ) == 9 );
}

TEST_CASE( "Phase 5A trace mesh enumerates boundary, patch-interface, and interior sides" )
{
    const MultiPatchSplineSpace h1 = makeTwoPatchSurface();

    const std::vector<TraceMeshInterface> boundaries = boundaryTraceMeshInterfaces( h1 );
    const std::vector<TraceMeshInterface> patch_interfaces = patchTraceMeshInterfaces( h1 );
    const std::vector<TraceMeshInterface> interiors = interiorTraceMeshInterfaces( h1 );
    const std::vector<TraceMeshInterface> broken_interiors =
        interiorTraceMeshInterfaces( h1, InteriorTraceMode::Broken );

    CHECK( boundaries.size() == 8 );
    CHECK( patch_interfaces.size() == 2 );
    CHECK( interiors.size() == 2 );
    CHECK( broken_interiors.empty() );

    for( const TraceMeshInterface& iface : boundaries )
    {
        CHECK( iface.type == TraceMeshEntityType::Boundary );
        CHECK( iface.isBoundary() );
        CHECK( not iface.permutation.has_value() );
        CHECK( traceInterfaceElement( h1, iface ).isBoundary() );
    }

    for( const TraceMeshInterface& iface : patch_interfaces )
    {
        CHECK( iface.type == TraceMeshEntityType::PatchInterface );
        REQUIRE_FALSE( iface.isBoundary() );
        REQUIRE( iface.permutation.has_value() );
        CHECK( iface.permutation.value() == TPPermutation::Flip1d );
        CHECK( iface.first.level == 0 );
        CHECK( iface.second->level == 0 );
        CHECK( iface.first.patch_id != iface.second->patch_id );
        CHECK( sharedTraceIdCount( traceInterfaceElement( h1, iface ) ) > 0 );
    }
    CHECK( aggregateSharedTraceIdCount( h1, patch_interfaces ) == 5 );

    for( const TraceMeshInterface& iface : interiors )
    {
        CHECK( iface.type == TraceMeshEntityType::Interior );
        REQUIRE_FALSE( iface.isBoundary() );
        CHECK( iface.first.patch_id == iface.second->patch_id );
        CHECK( sharedTraceIdCount( traceInterfaceElement( h1, iface ) ) == 3 );
    }
}

TEST_CASE( "Phase 5A broken interior filtering detects low-continuity reduced spaces" )
{
    const MultiPatchSplineSpace h1 = makeTwoPatchSurface();
    const VectorConformingMultiPatchSplineSpace hdiv = buildHDivMultiPatchSplineSpace( h1 );
    const MultiPatchSplineSpace l2 = buildL2MultiPatchSplineSpace( hdiv );

    const std::vector<TraceMeshInterface> interiors = interiorTraceMeshInterfaces( l2 );
    const std::vector<TraceMeshInterface> broken_interiors =
        interiorTraceMeshInterfaces( l2, InteriorTraceMode::Broken );

    REQUIRE( interiors.size() == 2 );
    CHECK( broken_interiors.size() == 2 );
    for( const TraceMeshInterface& iface : broken_interiors )
    {
        CHECK( sharedTraceIdCount( traceInterfaceElement( l2, iface ) ) == 0 );
    }
}

TEST_CASE( "Phase 5A trace mesh enumerates mixed-level hierarchical patch interfaces" )
{
    const HierarchicalMultiPatchSplineSpace h1 = makeMixedLevelTwoPatchHierarchy();
    const VectorConformingHierarchicalMultiPatchSplineSpace hdiv =
        buildHDivHierarchicalMultiPatchSplineSpace( h1 );

    const std::vector<TraceMeshInterface> boundaries = boundaryTraceMeshInterfaces( hdiv );
    const std::vector<TraceMeshInterface> patch_interfaces = patchTraceMeshInterfaces( hdiv );
    const std::vector<TraceMeshInterface> interiors = interiorTraceMeshInterfaces( hdiv );
    const std::vector<TraceMeshInterface> all =
        traceMeshInterfaces( hdiv, true, true, true );

    CHECK( boundaries.size() == 14 );
    REQUIRE( patch_interfaces.size() == 1 );
    CHECK( interiors.size() == 2 );
    CHECK( all.size() == boundaries.size() + patch_interfaces.size() + interiors.size() );

    const TraceMeshInterface& iface = patch_interfaces.front();
    CHECK( iface.type == TraceMeshEntityType::PatchInterface );
    REQUIRE_FALSE( iface.isBoundary() );
    CHECK( iface.first.level == 1 );
    CHECK( iface.second->level == 1 );
    CHECK( iface.first.patch_id == 0 );
    CHECK( iface.second->patch_id == 1 );
    REQUIRE( iface.permutation.has_value() );
    CHECK( iface.permutation.value() == TPPermutation::Flip1d );

    const TraceInterfaceElement trace = traceInterfaceElement( hdiv, iface );
    REQUIRE_FALSE( trace.isBoundary() );
    CHECK( sharedTraceIdCount( trace ) > 0 );
}
