#include <catch2/catch_test_macros.hpp>
#include <VectorConformingHierarchicalMultiPatchSplineSpace.hpp>
#include <MultiPatchSplineFactory.hpp>
#include <TraceExtraction.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParametricAtlas.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <IndexOperations.hpp>
#include <algorithm>
#include <set>

using namespace basis;
using namespace topology;

namespace
{
    using InternalConnectionsMap = MultiPatchCombinatorialMap::InternalConnectionsMap;
    using ConstituentSide = MultiPatchCombinatorialMap::ConstituentSide;
    using TPPermutation = MultiPatchCombinatorialMap::TPPermutation;

    std::shared_ptr<const TPSplineSpace> makePatch( const SmallVector<KnotVector, 3>& kvs,
                                                    const SmallVector<size_t, 3>& degrees )
    {
        return std::make_shared<const TPSplineSpace>( buildBSpline( kvs, degrees ) );
    }

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

    Cell patchTopCellAt( const MultiPatchCombinatorialMap& cmap,
                         const size_t patch_id,
                         const util::IndexVec& element_indices )
    {
        const Cell patch_cell = topCellAt( *cmap.constituents().at( patch_id ), element_indices );
        return Cell( cmap.toGlobalDart( patch_id, patch_cell.dart() ), cmap.dim() );
    }

    InternalConnectionsMap twoPatchConnection( const ElementSide& first_side,
                                               const ElementSide& second_side,
                                               const TPPermutation permutation )
    {
        const ConstituentSide first{ 0, first_side.sideId() };
        const ConstituentSide second{ 1, second_side.sideId() };
        return { { first, { permutation, second } }, { second, { permutation, first } } };
    }

    std::vector<Cell> topCells( const CombinatorialMap& cmap )
    {
        std::vector<Cell> out;
        iterateCellsWhile( cmap, cmap.dim(), [&]( const Cell& c ) {
            out.push_back( c );
            return true;
        } );
        return out;
    }

    std::vector<Cell> topCellsForPatch( const MultiPatchCombinatorialMap& cmap, const size_t patch_id )
    {
        std::vector<Cell> out;
        iterateCellsWhile( cmap, cmap.dim(), [&]( const Cell& c ) {
            const auto [cell_patch, local_d] = cmap.toLocalDart( c.dart() );
            if( cell_patch == patch_id ) out.emplace_back( cmap.toGlobalDart( patch_id, local_d ), c.dim() );
            return true;
        } );
        return out;
    }

    std::pair<HierarchicalMultiPatchSplineSpace, VectorConformingHierarchicalMultiPatchSplineSpace>
        makeAllFineTwoPatchHDivHierarchy()
    {
        const double ptol = 1e-10;
        const KnotVector coarse_kv( { 0, 0, 0, 1, 1, 1 }, ptol );
        const KnotVector fine_kv( { 0, 0, 0, 0.5, 1, 1, 1 }, ptol );

        const auto coarse_patch = makePatch( { coarse_kv, coarse_kv }, { 2, 2 } );
        const auto fine_patch = makePatch( { fine_kv, fine_kv }, { 2, 2 } );

        const ElementSide first_side( 0, false );
        const ElementSide second_side( 0, true );
        const InternalConnectionsMap connections =
            twoPatchConnection( first_side, second_side, TPPermutation::Flip1d );

        auto coarse_h1 = std::make_shared<const MultiPatchSplineSpace>(
            buildH1MultiPatchSplineSpace( { coarse_patch, coarse_patch }, connections ) );
        auto fine_h1 = std::make_shared<const MultiPatchSplineSpace>(
            buildH1MultiPatchSplineSpace( { fine_patch, fine_patch }, connections ) );

        const HierarchicalMultiPatchSplineSpace primal =
            buildHierarchicalSplineSpace(
                { coarse_h1, fine_h1 },
                { {}, topCells( fine_h1->basisComplex().parametricAtlas().cmap() ) } );

        return { primal, buildHDivHierarchicalMultiPatchSplineSpace( primal ) };
    }

    HierarchicalMultiPatchSplineSpace makeMixedLevelTwoPatchPrimalHierarchy()
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

    std::pair<HierarchicalMultiPatchSplineSpace, VectorConformingHierarchicalMultiPatchSplineSpace>
        makeAllFineTwoPatch3dHDivHierarchy()
    {
        const double ptol = 1e-10;
        const KnotVector coarse_kv( { 0, 0, 0, 1, 1, 1 }, ptol );
        const KnotVector fine_kv( { 0, 0, 0, 0.5, 1, 1, 1 }, ptol );

        const auto coarse_patch = makePatch( { coarse_kv, coarse_kv, coarse_kv }, { 2, 2, 2 } );
        const auto fine_patch = makePatch( { fine_kv, fine_kv, fine_kv }, { 2, 2, 2 } );

        const InternalConnectionsMap connections =
            twoPatchConnection( ElementSide( 0, false ), ElementSide( 0, true ), TPPermutation::ZeroToZero );

        auto coarse_h1 = std::make_shared<const MultiPatchSplineSpace>(
            buildH1MultiPatchSplineSpace( { coarse_patch, coarse_patch }, connections ) );
        auto fine_h1 = std::make_shared<const MultiPatchSplineSpace>(
            buildH1MultiPatchSplineSpace( { fine_patch, fine_patch }, connections ) );

        const HierarchicalMultiPatchSplineSpace primal =
            buildHierarchicalSplineSpace(
                { coarse_h1, fine_h1 },
                { {}, topCells( fine_h1->basisComplex().parametricAtlas().cmap() ) } );

        return { primal, buildHDivHierarchicalMultiPatchSplineSpace( primal ) };
    }

    std::vector<Cell> patchSideCells( const VectorConformingHierarchicalMultiPatchSplineSpace& ss,
                                      const size_t patch_id,
                                      const ElementSide& side )
    {
        const auto& cmap =
            dynamic_cast<const HierarchicalMultiPatchCombinatorialMap&>(
                ss.basisComplex().parametricAtlas().cmap() );
        const auto& constituent = *ss.constituents().at( patch_id );
        const auto& local_hier_cmap =
            dynamic_cast<const HierarchicalTPCombinatorialMap&>(
                constituent.basisComplex().parametricAtlas().cmap() );

        std::vector<Cell> out;
        iterateCellsWhile( cmap, cmap.dim(), [&]( const Cell& c ) {
            const auto [cell_patch, local_d] = cmap.dartRanges().toLocalDart( c.dart() );
            if( cell_patch != patch_id ) return true;

            const Cell local_c( local_d, c.dim() );
            const auto [level, level_d] = local_hier_cmap.unrefinedAncestorDartOfCell( local_c );
            const auto& level_cmap =
                dynamic_cast<const TPCombinatorialMap&>(
                    constituent.scalarBases().front()->refinementLevels().at( level )->basisComplex().parametricAtlas().cmap() );
            const FullyUnflattenedDart unflat = unflattenFull( level_cmap, level_d );
            const auto components = tensorProductComponentCMaps( level_cmap );
            const size_t axis_length =
                cellCount( *components.at( side.axis ), components.at( side.axis )->dim() );
            const size_t axis_index = unflat.unflat_darts.at( side.axis ).id();
            const bool on_side = side.lower ? axis_index == 0 : axis_index + 1 == axis_length;
            if( on_side ) out.push_back( c );
            return true;
        } );
        return out;
    }

    size_t localComponentIndex( const VectorConformingHierarchicalTPSplineSpace& ss, const FunctionId local_fid )
    {
        Eigen::Index offset = 0;
        for( size_t component = 0; component < ss.scalarBases().size(); component++ )
        {
            const Eigen::Index next_offset =
                offset + static_cast<Eigen::Index>( ss.scalarBases().at( component )->numFunctions() );
            if( local_fid.id() < next_offset ) return component;
            offset = next_offset;
        }
        throw std::runtime_error( "Local hierarchical vector function is outside component ranges." );
    }

    std::vector<std::set<Eigen::Index>> hierarchicalGlobalSideIdsByComponent(
        const VectorConformingHierarchicalMultiPatchSplineSpace& ss,
        const size_t patch_id,
        const ElementSide& side )
    {
        std::vector<std::set<Eigen::Index>> out( ss.numVectorComponents() );
        const auto& cmap =
            dynamic_cast<const HierarchicalMultiPatchCombinatorialMap&>(
                ss.basisComplex().parametricAtlas().cmap() );

        for( const Cell& global_cell : patchSideCells( ss, patch_id, side ) )
        {
            const auto [cell_patch, local_d] = cmap.dartRanges().toLocalDart( global_cell.dart() );
            REQUIRE( cell_patch == patch_id );
            const Cell local_cell( local_d, global_cell.dim() );
            const TraceSideData trace = traceSideData( *ss.constituents().at( patch_id ), local_cell, side );
            for( const FunctionId local_fid : trace.connectivity )
            {
                const size_t component = localComponentIndex( *ss.constituents().at( patch_id ), local_fid );
                out.at( component )
                    .insert( ss.functionIdMap().at( patch_id ).at( local_fid.id() ).first.id() );
            }
        }
        return out;
    }

    std::set<Eigen::Index> hierarchicalScalarGlobalSideIds( const HierarchicalMultiPatchSplineSpace& ss,
                                                            const size_t patch_id,
                                                            const ElementSide& side )
    {
        const auto& cmap =
            dynamic_cast<const HierarchicalMultiPatchCombinatorialMap&>(
                ss.basisComplex().parametricAtlas().cmap() );

        std::set<Eigen::Index> out;
        iterateCellsWhile( cmap, cmap.dim(), [&]( const Cell& c ) {
            const auto [cell_patch, local_d] = cmap.dartRanges().toLocalDart( c.dart() );
            if( cell_patch != patch_id ) return true;

            const Cell local_c( local_d, c.dim() );
            const auto& local_hier_cmap =
                dynamic_cast<const HierarchicalTPCombinatorialMap&>(
                    ss.constituents().at( patch_id )->basisComplex().parametricAtlas().cmap() );
            const auto [level, level_d] = local_hier_cmap.unrefinedAncestorDartOfCell( local_c );
            const auto& level_cmap =
                dynamic_cast<const TPCombinatorialMap&>(
                    ss.constituents().at( patch_id )
                        ->refinementLevels()
                        .at( level )
                        ->basisComplex()
                        .parametricAtlas()
                        .cmap() );
            const FullyUnflattenedDart unflat = unflattenFull( level_cmap, level_d );
            const auto components = tensorProductComponentCMaps( level_cmap );
            const size_t axis_length =
                cellCount( *components.at( side.axis ), components.at( side.axis )->dim() );
            const size_t axis_index = unflat.unflat_darts.at( side.axis ).id();
            const bool on_side = side.lower ? axis_index == 0 : axis_index + 1 == axis_length;
            if( not on_side ) return true;

            const TraceSideData trace = traceSideData( *ss.constituents().at( patch_id ), local_c, side );
            for( const FunctionId local_fid : trace.connectivity )
            {
                out.insert( ss.functionIdMap().at( patch_id ).at( local_fid.id() ).id() );
            }
            return true;
        } );
        return out;
    }

    size_t intersectionSize( const std::set<Eigen::Index>& first, const std::set<Eigen::Index>& second )
    {
        std::set<Eigen::Index> out;
        std::set_intersection(
            first.begin(), first.end(), second.begin(), second.end(), std::inserter( out, out.end() ) );
        return out.size();
    }

    Eigen::Index rowOf( const std::vector<FunctionId>& conn, const FunctionId fid )
    {
        for( Eigen::Index i = 0; i < static_cast<Eigen::Index>( conn.size() ); i++ )
        {
            if( conn.at( static_cast<size_t>( i ) ).id() == fid.id() ) return i;
        }
        throw std::runtime_error( "Function id was not found in trace connectivity." );
    }

    std::pair<size_t, TraceSideData> levelTraceSideData(
        const VectorConformingHierarchicalMultiPatchSplineSpace& ss,
        const Cell& hierarchical_cell,
        const ElementSide& side )
    {
        const auto& cmap =
            dynamic_cast<const HierarchicalMultiPatchCombinatorialMap&>(
                ss.basisComplex().parametricAtlas().cmap() );
        const auto [level, level_d] = cmap.unrefinedAncestorDartOfCell( hierarchical_cell );
        return { level,
                 traceSideData(
                     *ss.refinementLevels().at( level ), Cell( level_d, hierarchical_cell.dim() ), side ) };
    }

    param::ParentPoint centerPoint( const basis::SplineSpace& ss, const Cell& c )
    {
        const size_t dim = ss.basisComplex().parametricAtlas().cmap().dim();
        Vector3dMax pt( dim );
        pt.setConstant( 0.5 );
        return param::ParentPoint(
            ss.basisComplex().parametricAtlas().parentDomain( c ),
            pt,
            param::BaryCoordIsZeroVec( 2 * dim, false ) );
    }

    void checkEvaluatorOnAllLeafElements( const basis::SplineSpace& ss )
    {
        eval::SplineSpaceEvaluator evaler( ss, 0 );
        iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), ss.basisComplex().parametricAtlas().cmap().dim(), [&]( const Cell& c ) {
            evaler.localizeElement( c );
            evaler.localizePoint( centerPoint( ss, c ) );
            const Eigen::MatrixXd basis = evaler.evaluateBasis();
            CHECK( basis.rows() == static_cast<Eigen::Index>( ss.connectivity( c ).size() ) );
            CHECK( basis.cols() == static_cast<Eigen::Index>( ss.numVectorComponents() ) );
            CHECK( basis.array().isFinite().all() );
            return true;
        } );
    }
}

TEST_CASE( "Phase 4a builds hierarchical H(div) multipatch space with all fine leaves" )
{
    const auto [primal, hdiv] = makeAllFineTwoPatchHDivHierarchy();

    REQUIRE( hdiv.refinementLevels().size() == 2 );
    CHECK( hdiv.numFunctions() == hdiv.refinementLevels().at( 1 )->numFunctions() );
    CHECK( hdiv.numVectorComponents() == 2 );

    const std::vector<std::pair<size_t, FunctionId>> global_level_map =
        vectorFidsToRefinementLevelFids( hdiv );
    REQUIRE( global_level_map.size() == hdiv.numFunctions() );
    for( const auto& [level, fid] : global_level_map )
    {
        CHECK( level == 1 );
        CHECK( fid.id() >= 0 );
    }

    for( size_t patch_id = 0; patch_id < hdiv.constituents().size(); patch_id++ )
    {
        const auto patch_level_map = patchVectorFidsToRefinementLevelFids( hdiv, patch_id );
        CHECK( patch_level_map.size() == hdiv.functionIdMap().at( patch_id ).size() );
        CHECK( patch_level_map.size() == hdiv.constituents().at( patch_id )->numFunctions() );
        for( const auto& entry : patch_level_map )
        {
            CHECK( entry.level == 1 );
        }
    }

    (void)primal;
}

TEST_CASE( "Phase 4a hierarchical H(div) shares only normal component across patch interfaces" )
{
    const auto [primal, hdiv] = makeAllFineTwoPatchHDivHierarchy();
    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );

    const auto first_ids = hierarchicalGlobalSideIdsByComponent( hdiv, 0, first_side );
    const auto second_ids = hierarchicalGlobalSideIdsByComponent( hdiv, 1, second_side );
    REQUIRE( first_ids.size() == 2 );
    REQUIRE( second_ids.size() == 2 );

    CHECK( intersectionSize( first_ids.at( 0 ), second_ids.at( 0 ) ) > 0 );
    CHECK( intersectionSize( first_ids.at( 1 ), second_ids.at( 1 ) ) == 0 );

    (void)primal;
}

TEST_CASE( "Phase 4a hierarchical vector extraction applies interface orientations to trace rows" )
{
    const auto [primal, hdiv] = makeAllFineTwoPatchHDivHierarchy();
    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );

    const std::vector<Cell> first_cells = patchSideCells( hdiv, 0, first_side );
    const std::vector<Cell> second_cells = patchSideCells( hdiv, 1, second_side );
    REQUIRE_FALSE( first_cells.empty() );
    REQUIRE( first_cells.size() == second_cells.size() );

    const Cell first_cell = first_cells.front();
    const Cell second_cell = second_cells.back();
    const TraceSideData first_trace = traceSideData( hdiv, first_cell, first_side );
    const TraceSideData second_trace = traceSideData( hdiv, second_cell, second_side );
    const auto [first_level, first_level_trace] = levelTraceSideData( hdiv, first_cell, first_side );
    const auto [second_level, second_level_trace] = levelTraceSideData( hdiv, second_cell, second_side );
    const std::vector<std::pair<size_t, FunctionId>> global_level_map =
        vectorFidsToRefinementLevelFids( hdiv );

    bool checked_shared_row = false;
    for( const FunctionId fid : first_trace.connectivity )
    {
        if( std::ranges::none_of( second_trace.connectivity, [&]( const FunctionId other ) {
                return other.id() == fid.id();
            } ) )
            continue;

        const Eigen::Index first_row = rowOf( first_trace.connectivity, fid );
        const Eigen::Index second_row = rowOf( second_trace.connectivity, fid );
        const auto [level_from_map, level_fid] = global_level_map.at( fid.id() );
        REQUIRE( level_from_map == first_level );
        REQUIRE( level_from_map == second_level );
        const Eigen::Index first_level_row = rowOf( first_level_trace.connectivity, level_fid );
        const Eigen::Index second_level_row = rowOf( second_level_trace.connectivity, level_fid );

        CHECK( ( first_trace.extraction.row( first_row ) -
                 first_level_trace.extraction.row( first_level_row ) )
                   .norm() < 1e-10 );
        CHECK( ( second_trace.extraction.row( second_row ) -
                 second_level_trace.extraction.row( second_level_row ) )
                   .norm() < 1e-10 );
        checked_shared_row = true;
    }

    CHECK( checked_shared_row );
    (void)primal;
}

TEST_CASE( "Phase 4a hierarchical H(curl) can be built over the same primal hierarchy" )
{
    const auto [primal, hdiv] = makeAllFineTwoPatchHDivHierarchy();
    const VectorConformingHierarchicalMultiPatchSplineSpace hcurl =
        buildHCurlHierarchicalMultiPatchSplineSpace( primal );

    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );
    const auto first_ids = hierarchicalGlobalSideIdsByComponent( hcurl, 0, first_side );
    const auto second_ids = hierarchicalGlobalSideIdsByComponent( hcurl, 1, second_side );

    REQUIRE( first_ids.size() == 2 );
    CHECK( intersectionSize( first_ids.at( 0 ), second_ids.at( 0 ) ) == 0 );
    CHECK( intersectionSize( first_ids.at( 1 ), second_ids.at( 1 ) ) > 0 );

    const auto& cmap = hcurl.basisComplex().parametricAtlas().cmap();
    const ParentBasis pb = hcurl.basisComplex().parentBasis( topCells( cmap ).front() );
    CHECK( pb.mBasisGroups.size() == 1 );

    (void)hdiv;
}

TEST_CASE( "Phase 4a mixed-level 2D hierarchy keeps coarse away-from-interface functions active" )
{
    const HierarchicalMultiPatchSplineSpace primal = makeMixedLevelTwoPatchPrimalHierarchy();
    const VectorConformingHierarchicalMultiPatchSplineSpace hdiv =
        buildHDivHierarchicalMultiPatchSplineSpace( primal );
    const VectorConformingHierarchicalMultiPatchSplineSpace hcurl =
        buildHCurlHierarchicalMultiPatchSplineSpace( primal );

    CHECK( hdiv.numFunctions() != hdiv.refinementLevels().at( 1 )->numFunctions() );
    CHECK( hcurl.numFunctions() != hcurl.refinementLevels().at( 1 )->numFunctions() );

    const std::vector<std::pair<size_t, FunctionId>> hdiv_level_map =
        vectorFidsToRefinementLevelFids( hdiv );
    const bool hdiv_has_level_0 =
        std::ranges::any_of( hdiv_level_map, []( const auto& pr ) { return pr.first == 0; } );
    const bool hdiv_has_level_1 =
        std::ranges::any_of( hdiv_level_map, []( const auto& pr ) { return pr.first == 1; } );
    CHECK( hdiv_has_level_0 );
    CHECK( hdiv_has_level_1 );

    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );

    const auto hdiv_first = hierarchicalGlobalSideIdsByComponent( hdiv, 0, first_side );
    const auto hdiv_second = hierarchicalGlobalSideIdsByComponent( hdiv, 1, second_side );
    REQUIRE( hdiv_first.size() == 2 );
    CHECK( intersectionSize( hdiv_first.at( 0 ), hdiv_second.at( 0 ) ) > 0 );
    CHECK( intersectionSize( hdiv_first.at( 1 ), hdiv_second.at( 1 ) ) == 0 );

    const auto hcurl_first = hierarchicalGlobalSideIdsByComponent( hcurl, 0, first_side );
    const auto hcurl_second = hierarchicalGlobalSideIdsByComponent( hcurl, 1, second_side );
    REQUIRE( hcurl_first.size() == 2 );
    CHECK( intersectionSize( hcurl_first.at( 0 ), hcurl_second.at( 0 ) ) == 0 );
    CHECK( intersectionSize( hcurl_first.at( 1 ), hcurl_second.at( 1 ) ) > 0 );
}

TEST_CASE( "Phase 4a rejects nonmatching hierarchical patch interfaces" )
{
    const double ptol = 1e-10;
    const KnotVector coarse_kv( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector fine_kv( { 0, 0, 0, 0.5, 1, 1, 1 }, ptol );

    const auto coarse_patch = makePatch( { coarse_kv, coarse_kv }, { 2, 2 } );
    const auto fine_patch = makePatch( { fine_kv, fine_kv }, { 2, 2 } );
    const InternalConnectionsMap connections =
        twoPatchConnection( ElementSide( 0, false ), ElementSide( 0, true ), TPPermutation::Flip1d );

    auto coarse_h1 = std::make_shared<const MultiPatchSplineSpace>(
        buildH1MultiPatchSplineSpace( { coarse_patch, coarse_patch }, connections ) );
    auto fine_h1 = std::make_shared<const MultiPatchSplineSpace>(
        buildH1MultiPatchSplineSpace( { fine_patch, fine_patch }, connections ) );

    auto coarse_hdiv = std::make_shared<const VectorConformingMultiPatchSplineSpace>(
        buildHDivMultiPatchSplineSpace( *coarse_h1 ) );
    auto fine_hdiv = std::make_shared<const VectorConformingMultiPatchSplineSpace>(
        buildHDivMultiPatchSplineSpace( *fine_h1 ) );

    const auto& coarse_cmap =
        dynamic_cast<const MultiPatchCombinatorialMap&>(
            coarse_h1->basisComplex().parametricAtlas().cmap() );
    const auto& fine_cmap =
        dynamic_cast<const MultiPatchCombinatorialMap&>(
            fine_h1->basisComplex().parametricAtlas().cmap() );

    CHECK_THROWS_AS(
        buildHierarchicalSplineSpace(
            { coarse_hdiv, fine_hdiv },
            { topCellsForPatch( coarse_cmap, 1 ), topCellsForPatch( fine_cmap, 0 ) } ),
        std::invalid_argument );
}

TEST_CASE( "Phase 4a all-fine 3D hierarchy has exact H(div) and H(curl) face sharing" )
{
    const auto [primal, hdiv] = makeAllFineTwoPatch3dHDivHierarchy();
    const VectorConformingHierarchicalMultiPatchSplineSpace hcurl =
        buildHCurlHierarchicalMultiPatchSplineSpace( primal );

    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );

    const auto hdiv_first = hierarchicalGlobalSideIdsByComponent( hdiv, 0, first_side );
    const auto hdiv_second = hierarchicalGlobalSideIdsByComponent( hdiv, 1, second_side );
    REQUIRE( hdiv_first.size() == 3 );
    CHECK( intersectionSize( hdiv_first.at( 0 ), hdiv_second.at( 0 ) ) == 9 );
    CHECK( intersectionSize( hdiv_first.at( 1 ), hdiv_second.at( 1 ) ) == 0 );
    CHECK( intersectionSize( hdiv_first.at( 2 ), hdiv_second.at( 2 ) ) == 0 );

    const auto hcurl_first = hierarchicalGlobalSideIdsByComponent( hcurl, 0, first_side );
    const auto hcurl_second = hierarchicalGlobalSideIdsByComponent( hcurl, 1, second_side );
    REQUIRE( hcurl_first.size() == 3 );
    CHECK( intersectionSize( hcurl_first.at( 0 ), hcurl_second.at( 0 ) ) == 0 );
    CHECK( intersectionSize( hcurl_first.at( 1 ), hcurl_second.at( 1 ) ) == 12 );
    CHECK( intersectionSize( hcurl_first.at( 2 ), hcurl_second.at( 2 ) ) == 12 );

    const std::vector<std::pair<size_t, FunctionId>> global_level_map =
        vectorFidsToRefinementLevelFids( hdiv );
    REQUIRE( global_level_map.size() == hdiv.numFunctions() );
    for( size_t patch_id = 0; patch_id < hdiv.constituents().size(); patch_id++ )
    {
        const auto patch_level_map = patchVectorFidsToRefinementLevelFids( hdiv, patch_id );
        CHECK( patch_level_map.size() == hdiv.functionIdMap().at( patch_id ).size() );
        CHECK( patch_level_map.size() == hdiv.constituents().at( patch_id )->numFunctions() );
    }
}

TEST_CASE( "Phase 4b builds discontinuous hierarchical L2 from hierarchical H(div)" )
{
    const auto [primal, hdiv] = makeAllFineTwoPatchHDivHierarchy();
    const HierarchicalMultiPatchSplineSpace l2 = buildL2HierarchicalMultiPatchSplineSpace( hdiv );

    const MultiPatchSplineSpace level_l2 = buildL2MultiPatchSplineSpace( *hdiv.refinementLevels().at( 1 ) );
    CHECK( l2.numFunctions() == level_l2.numFunctions() );

    const std::vector<std::pair<size_t, FunctionId>> level_map = fidsToRefinementLevelFids( l2 );
    REQUIRE( level_map.size() == l2.numFunctions() );
    for( const auto& [level, fid] : level_map )
    {
        CHECK( level == 1 );
        CHECK( fid.id() >= 0 );
    }

    const ElementSide first_side( 0, false );
    const ElementSide second_side( 0, true );
    CHECK( intersectionSize( hierarchicalScalarGlobalSideIds( l2, 0, first_side ),
                             hierarchicalScalarGlobalSideIds( l2, 1, second_side ) ) == 0 );

    (void)primal;
}

TEST_CASE( "Phase 4a hierarchical vector spaces evaluate on all leaf element centers" )
{
    const HierarchicalMultiPatchSplineSpace mixed_primal = makeMixedLevelTwoPatchPrimalHierarchy();
    const VectorConformingHierarchicalMultiPatchSplineSpace mixed_hdiv =
        buildHDivHierarchicalMultiPatchSplineSpace( mixed_primal );
    const VectorConformingHierarchicalMultiPatchSplineSpace mixed_hcurl =
        buildHCurlHierarchicalMultiPatchSplineSpace( mixed_primal );
    checkEvaluatorOnAllLeafElements( mixed_hdiv );
    checkEvaluatorOnAllLeafElements( mixed_hcurl );

    const auto [primal_3d, hdiv_3d] = makeAllFineTwoPatch3dHDivHierarchy();
    const VectorConformingHierarchicalMultiPatchSplineSpace hcurl_3d =
        buildHCurlHierarchicalMultiPatchSplineSpace( primal_3d );
    checkEvaluatorOnAllLeafElements( hdiv_3d );
    checkEvaluatorOnAllLeafElements( hcurl_3d );
}
