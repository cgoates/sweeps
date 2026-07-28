#include <VectorConformingHierarchicalMultiPatchSplineSpace.hpp>
#include <MultiPatchSplineFactory.hpp>
#include <CombinatorialMapMethods.hpp>
#include <TraceExtraction.hpp>
#include <TraceMesh.hpp>
#include <HierarchicalMultiPatchCombinatorialMap.hpp>
#include <MultiPatchBasisComplex.hpp>
#include <HierarchicalMultiPatchParametricAtlas.hpp>
#include <IndexOperations.hpp>
#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

using namespace basis;

namespace
{
    using SideIndex = std::vector<size_t>;

    struct VectorActiveFuncs
    {
        std::vector<std::vector<std::vector<FunctionId>>> active_vec_funcs;
        std::vector<std::vector<std::pair<FunctionId, bool>>> patch_to_global;
        size_t num_active_funcs;
    };

    std::vector<size_t> componentOffsets( const VectorConformingTPSplineSpace& ss )
    {
        std::vector<size_t> offsets;
        offsets.reserve( ss.scalarTPBases().size() + 1 );
        offsets.push_back( 0 );
        for( const auto& scalar_basis : ss.scalarTPBases() )
        {
            offsets.push_back( offsets.back() + scalar_basis->numFunctions() );
        }
        if( offsets.back() != ss.numFunctions() )
            throw std::runtime_error( "Vector component offsets do not span the full space." );
        return offsets;
    }

    bool inComponentRange( const FunctionId fid, const std::vector<size_t>& offsets, const size_t component )
    {
        return offsets.at( component ) <= static_cast<size_t>( fid.id() ) and
               static_cast<size_t>( fid.id() ) < offsets.at( component + 1 );
    }

    SideIndex sideLengths( const topology::TPCombinatorialMap& cmap, const ElementSide& side )
    {
        const SmallVector<std::shared_ptr<const topology::CombinatorialMap1d>, 3> components =
            topology::tensorProductComponentCMaps( cmap );

        SideIndex out;
        for( size_t axis = 0; axis < cmap.dim(); axis++ )
        {
            if( axis == side.axis ) continue;
            const auto& component = components.at( axis );
            out.push_back( topology::cellCount( *component, component->dim() ) );
        }
        return out;
    }

    std::set<SideIndex> sideLeafIndices( const topology::MultiPatchCombinatorialMap& level_cmap,
                                         const std::vector<topology::Cell>& leaf_elements,
                                         const size_t patch_id,
                                         const ElementSide& side )
    {
        const auto& patch_cmap = *level_cmap.constituents().at( patch_id );
        const SmallVector<std::shared_ptr<const topology::CombinatorialMap1d>, 3> components =
            topology::tensorProductComponentCMaps( patch_cmap );
        const size_t normal_length =
            topology::cellCount( *components.at( side.axis ), components.at( side.axis )->dim() );

        std::set<SideIndex> out;
        for( const topology::Cell& leaf : leaf_elements )
        {
            if( leaf.dim() != level_cmap.dim() ) continue;
            const auto [leaf_patch, local_d] = level_cmap.toLocalDart( leaf.dart() );
            if( leaf_patch != patch_id ) continue;

            const topology::FullyUnflattenedDart unflat = topology::unflattenFull( patch_cmap, local_d );
            const size_t normal_index = unflat.unflat_darts.at( side.axis ).id();
            const bool on_side = side.lower ? normal_index == 0 : normal_index + 1 == normal_length;
            if( not on_side ) continue;

            SideIndex side_index;
            for( size_t axis = 0; axis < patch_cmap.dim(); axis++ )
            {
                if( axis != side.axis ) side_index.push_back( unflat.unflat_darts.at( axis ).id() );
            }
            out.insert( side_index );
        }
        return out;
    }

    void validateMatchingInterfaces(
        const topology::HierarchicalMultiPatchCombinatorialMap& cmap,
        const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>>& refinement_levels )
    {
        const std::vector<std::vector<topology::Cell>> leaves = topology::leafElements( cmap );
        if( leaves.size() != refinement_levels.size() )
            throw std::runtime_error( "Hierarchical multipatch leaf data does not match refinement levels." );

        for( size_t level_ii = 0; level_ii < refinement_levels.size(); level_ii++ )
        {
            const topology::MultiPatchCombinatorialMap& level_cmap = *cmap.refinementLevels().at( level_ii );
            for( const auto& [first, connection] : level_cmap.connections() )
            {
                const topology::MultiPatchCombinatorialMap::ConstituentSide& second = connection.second;
                if( not( first < second ) ) continue;

                const ElementSide first_side = elementSideFromId( first.side_id );
                const ElementSide second_side = elementSideFromId( second.side_id );
                const SideIndex first_lengths =
                    sideLengths( *level_cmap.constituents().at( first.constituent_id ), first_side );

                const std::set<SideIndex> first_indices =
                    sideLeafIndices( level_cmap, leaves.at( level_ii ), first.constituent_id, first_side );
                const std::set<SideIndex> second_indices =
                    sideLeafIndices( level_cmap, leaves.at( level_ii ), second.constituent_id, second_side );

                std::set<SideIndex> permuted_first_indices;
                for( const SideIndex& index : first_indices )
                {
                    permuted_first_indices.insert( permuteTraceSideIndex( index, first_lengths, connection.first ) );
                }

                if( permuted_first_indices != second_indices )
                    throw std::invalid_argument(
                        "VectorConformingHierarchicalMultiPatchSplineSpace requires matching leaf elements across patch interfaces." );
            }
        }
    }

    VectorActiveFuncs vectorActiveFuncs(
        const topology::HierarchicalMultiPatchCombinatorialMap& cmap,
        const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>>& refinement_levels )
    {
        const std::vector<std::vector<topology::Cell>> leaf_elements = topology::leafElements( cmap );
        std::vector<std::vector<std::vector<FunctionId>>> active_funcs(
            refinement_levels.front()->subSpaces().size(),
            std::vector<std::vector<FunctionId>>( refinement_levels.size() ) );
        std::vector<std::map<FunctionId, FunctionId>> level_to_hierarchical_fids( refinement_levels.size() );

        size_t n_active_funcs = 0;
        for( size_t level_ii = 0; level_ii < refinement_levels.size(); level_ii++ )
        {
            const auto& level_ss = refinement_levels.at( level_ii );
            const auto& level_cmap = *cmap.refinementLevels().at( level_ii );
            std::set<topology::Cell> cells_to_prune;
            std::set<FunctionId> level_active_funcs;
            for( const topology::Cell& leaf_elem : leaf_elements.at( level_ii ) )
            {
                const std::vector<FunctionId> conn = level_ss->connectivity( leaf_elem );
                level_active_funcs.insert( conn.begin(), conn.end() );
                topology::iterateAdjacentCells( level_cmap, leaf_elem, cmap.dim() - 1, [&]( const topology::Cell& c ) {
                    const auto maybe_phi = topology::phi( level_cmap, cmap.dim(), c.dart() );
                    if( maybe_phi )
                    {
                        bool ancestor_leaf = false;
                        topology::iterateDartsOfCell(
                            level_cmap,
                            topology::Cell( maybe_phi.value(), cmap.dim() ),
                            [&]( const topology::Dart& d ) {
                                cmap.iterateAncestors( cmap.toGlobalDart( level_ii, d ), [&]( const topology::Dart& ancestor ) {
                                    if( cmap.isUnrefinedLeafDart( ancestor ) )
                                    {
                                        ancestor_leaf = true;
                                    }
                                    return not ancestor_leaf;
                                } );
                                return not ancestor_leaf;
                            } );
                        if( ancestor_leaf )
                        {
                            cells_to_prune.insert( topology::Cell( maybe_phi.value(), cmap.dim() ) );
                        }
                    }
                    return true;
                } );
            }

            for( const topology::Cell& prune_cell : cells_to_prune )
            {
                const std::vector<FunctionId> conn = level_ss->connectivity( prune_cell );
                for( const FunctionId& fid : conn ) level_active_funcs.erase( fid );
            }

            size_t level_active_index = 0;
            for( const FunctionId& fid : level_active_funcs )
            {
                level_to_hierarchical_fids.at( level_ii ).emplace(
                    fid, FunctionId( n_active_funcs + level_active_index++ ) );
            }

            for( size_t patch_ii = 0; patch_ii < level_ss->subSpaces().size(); patch_ii++ )
            {
                auto& patch_level_active_funcs = active_funcs.at( patch_ii ).at( level_ii );
                const auto& patch_funcs_to_mp_funcs = level_ss->functionIdMap().at( patch_ii );
                for( size_t patch_fid = 0; patch_fid < patch_funcs_to_mp_funcs.size(); patch_fid++ )
                {
                    if( level_active_funcs.contains( patch_funcs_to_mp_funcs.at( patch_fid ).first ) )
                    {
                        patch_level_active_funcs.push_back( FunctionId( patch_fid ) );
                    }
                }
            }
            n_active_funcs += level_active_funcs.size();
        }

        std::vector<std::vector<std::pair<FunctionId, bool>>> patch_to_global(
            refinement_levels.front()->subSpaces().size() );
        for( size_t patch_ii = 0; patch_ii < patch_to_global.size(); patch_ii++ )
        {
            for( size_t component = 0; component < refinement_levels.front()->numVectorComponents(); component++ )
            {
                for( size_t level_ii = 0; level_ii < refinement_levels.size(); level_ii++ )
                {
                    const auto& patch_level = *refinement_levels.at( level_ii )->subSpaces().at( patch_ii );
                    const std::vector<size_t> offsets = componentOffsets( patch_level );
                    const auto& patch_level_map = refinement_levels.at( level_ii )->functionIdMap().at( patch_ii );
                    for( const FunctionId& patch_fid : active_funcs.at( patch_ii ).at( level_ii ) )
                    {
                        if( not inComponentRange( patch_fid, offsets, component ) ) continue;

                        const auto& [level_global_fid, orientation] = patch_level_map.at( patch_fid.id() );
                        const auto it = level_to_hierarchical_fids.at( level_ii ).find( level_global_fid );
                        if( it == level_to_hierarchical_fids.at( level_ii ).end() )
                            throw std::runtime_error( "Active patch vector function was not active globally." );

                        patch_to_global.at( patch_ii ).push_back( { it->second, orientation } );
                    }
                }
            }
        }

        return { active_funcs, patch_to_global, n_active_funcs };
    }

    std::vector<std::shared_ptr<const VectorConformingHierarchicalTPSplineSpace>>
        initializeConstituents(
            const VectorConformingBasisComplex& bc,
            const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>>& refinement_levels,
            const std::vector<std::vector<std::vector<FunctionId>>>& active_funcs )
    {
        const auto primal_hier_bc =
            std::dynamic_pointer_cast<const HierarchicalMultiPatchBasisComplex>( bc.primalComplexPtr() );
        if( primal_hier_bc == nullptr )
            throw std::invalid_argument(
                "VectorConformingHierarchicalMultiPatchSplineSpace requires a HierarchicalMultiPatchBasisComplex." );

        if( primal_hier_bc->refinementLevels().size() != refinement_levels.size() )
            throw std::invalid_argument( "Inconsistent refinement levels and basis complex." );

        std::vector<std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>> constituent_levels(
            primal_hier_bc->constituents().size() );
        for( const auto& level : refinement_levels )
        {
            if( level->subSpaces().size() != primal_hier_bc->constituents().size() )
                throw std::invalid_argument( "All vector multipatch levels must have the same number of patches." );
            for( size_t patch_ii = 0; patch_ii < primal_hier_bc->constituents().size(); patch_ii++ )
            {
                constituent_levels.at( patch_ii ).push_back( level->subSpaces().at( patch_ii ) );
            }
        }

        std::vector<std::shared_ptr<const VectorConformingHierarchicalTPSplineSpace>> out;
        out.reserve( primal_hier_bc->constituents().size() );
        for( size_t patch_ii = 0; patch_ii < primal_hier_bc->constituents().size(); patch_ii++ )
        {
            const auto patch_bc =
                std::make_shared<const VectorConformingBasisComplex>(
                    primal_hier_bc->constituents().at( patch_ii ), bc.conformingType() );
            out.push_back( std::make_shared<const VectorConformingHierarchicalTPSplineSpace>(
                patch_bc, constituent_levels.at( patch_ii ), active_funcs.at( patch_ii ) ) );
        }
        return out;
    }
}

namespace basis
{
    VectorConformingHierarchicalMultiPatchSplineSpace::VectorConformingHierarchicalMultiPatchSplineSpace(
        const std::shared_ptr<const VectorConformingBasisComplex>& bc,
        const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>>& refinement_levels )
        : mBasisComplex( bc ), mRefinementLevels( refinement_levels )
    {
        if( refinement_levels.empty() )
            throw std::invalid_argument( "VectorConformingHierarchicalMultiPatchSplineSpace requires refinement levels." );

        const auto primal_hier_bc =
            std::dynamic_pointer_cast<const HierarchicalMultiPatchBasisComplex>( bc->primalComplexPtr() );
        if( primal_hier_bc == nullptr )
            throw std::invalid_argument(
                "VectorConformingHierarchicalMultiPatchSplineSpace requires a HierarchicalMultiPatchBasisComplex." );

        validateMatchingInterfaces( primal_hier_bc->parametricAtlas().cmap(), refinement_levels );

        const VectorActiveFuncs active_funcs =
            vectorActiveFuncs( primal_hier_bc->parametricAtlas().cmap(), refinement_levels );
        mConstituents = initializeConstituents( *bc, refinement_levels, active_funcs.active_vec_funcs );
        mFuncIds = active_funcs.patch_to_global;
        mNumActiveFuncs = active_funcs.num_active_funcs;

        if( mFuncIds.size() != mConstituents.size() )
            throw std::runtime_error( "Patch-to-global vector function map has inconsistent patch count." );
        for( size_t patch_ii = 0; patch_ii < mFuncIds.size(); patch_ii++ )
        {
            if( mFuncIds.at( patch_ii ).size() != mConstituents.at( patch_ii )->numFunctions() )
                throw std::runtime_error( "Patch-to-global vector function map has inconsistent local function count." );
        }
    }

    const VectorConformingBasisComplex& VectorConformingHierarchicalMultiPatchSplineSpace::basisComplex() const
    {
        return *mBasisComplex;
    }

    Eigen::MatrixXd
        VectorConformingHierarchicalMultiPatchSplineSpace::extractionOperator( const topology::Cell& c ) const
    {
        const auto& cmap =
            dynamic_cast<const topology::HierarchicalMultiPatchCombinatorialMap&>(
                mBasisComplex->parametricAtlas().cmap() );
        const auto [patch_ii, constituent_d] = cmap.dartRanges().toLocalDart( c.dart() );
        const topology::Cell constituent_c( constituent_d, c.dim() );
        const std::vector<FunctionId> local_conn = mConstituents.at( patch_ii )->connectivity( constituent_c );

        Eigen::VectorXd orientations( local_conn.size() );
        for( size_t i = 0; i < local_conn.size(); i++ )
        {
            const auto& [global_fid, orientation] = mFuncIds.at( patch_ii ).at( local_conn.at( i ).id() );
            (void)global_fid;
            orientations( static_cast<Eigen::Index>( i ) ) = orientation ? 1.0 : -1.0;
        }

        return orientations.asDiagonal() * mConstituents.at( patch_ii )->extractionOperator( constituent_c );
    }

    std::vector<FunctionId>
        VectorConformingHierarchicalMultiPatchSplineSpace::connectivity( const topology::Cell& c ) const
    {
        const auto& cmap =
            dynamic_cast<const topology::HierarchicalMultiPatchCombinatorialMap&>(
                mBasisComplex->parametricAtlas().cmap() );
        const auto [patch_ii, constituent_d] = cmap.dartRanges().toLocalDart( c.dart() );
        const topology::Cell constituent_c( constituent_d, c.dim() );
        const std::vector<FunctionId> local_conn = mConstituents.at( patch_ii )->connectivity( constituent_c );

        std::vector<FunctionId> global_conn;
        global_conn.reserve( local_conn.size() );
        std::transform( local_conn.begin(),
                        local_conn.end(),
                        std::back_inserter( global_conn ),
                        [&]( const FunctionId& local_fid ) {
                            return mFuncIds.at( patch_ii ).at( local_fid.id() ).first;
                        } );

        return global_conn;
    }

    size_t VectorConformingHierarchicalMultiPatchSplineSpace::numFunctions() const
    {
        return mNumActiveFuncs;
    }

    size_t VectorConformingHierarchicalMultiPatchSplineSpace::numVectorComponents() const
    {
        return mConstituents.front()->numVectorComponents();
    }

    VectorConformingHierarchicalMultiPatchSplineSpace
        buildHierarchicalSplineSpace(
            const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>>& refinement_levels,
            const std::vector<std::vector<topology::Cell>>& leaf_elements )
    {
        if( refinement_levels.empty() )
            throw std::invalid_argument( "Cannot build hierarchical vector multipatch space with no levels." );

        std::vector<std::shared_ptr<const MultiPatchBasisComplex>> bc_levels;
        bc_levels.reserve( refinement_levels.size() );
        std::vector<std::shared_ptr<const param::MultiPatchParametricAtlas>> atlas_levels;
        atlas_levels.reserve( refinement_levels.size() );
        std::vector<std::shared_ptr<const topology::MultiPatchCombinatorialMap>> cmap_levels;
        cmap_levels.reserve( refinement_levels.size() );

        const ConformingType conforming_type = refinement_levels.front()->basisComplex().conformingType();
        for( const auto& ss : refinement_levels )
        {
            if( ss->basisComplex().conformingType() != conforming_type )
                throw std::invalid_argument( "All vector multipatch levels must use the same conforming type." );

            const auto primal_mp_bc =
                std::dynamic_pointer_cast<const MultiPatchBasisComplex>( ss->basisComplex().primalComplexPtr() );
            if( primal_mp_bc == nullptr )
                throw std::invalid_argument( "Vector multipatch refinement level does not wrap a MultiPatchBasisComplex." );

            bc_levels.push_back( primal_mp_bc );
            atlas_levels.push_back( primal_mp_bc->parametricAtlasPtr() );
            cmap_levels.push_back( atlas_levels.back()->cmapPtr() );
        }

        const auto cmap =
            std::make_shared<const topology::HierarchicalMultiPatchCombinatorialMap>( cmap_levels, leaf_elements );
        const auto atlas =
            std::make_shared<const param::HierarchicalMultiPatchParametricAtlas>( cmap, atlas_levels );
        const auto primal_bc =
            std::make_shared<const HierarchicalMultiPatchBasisComplex>( atlas, bc_levels );
        const auto vector_bc =
            std::make_shared<const VectorConformingBasisComplex>( primal_bc, conforming_type );

        return VectorConformingHierarchicalMultiPatchSplineSpace( vector_bc, refinement_levels );
    }

    VectorConformingHierarchicalMultiPatchSplineSpace
        buildVectorConformingHierarchicalMultiPatchSplineSpace( const HierarchicalMultiPatchSplineSpace& primal,
                                                                const ConformingType conforming_type )
    {
        std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>> vector_levels;
        vector_levels.reserve( primal.refinementLevels().size() );
        for( const auto& level : primal.refinementLevels() )
        {
            vector_levels.push_back(
                std::make_shared<const VectorConformingMultiPatchSplineSpace>(
                    buildVectorConformingMultiPatchSplineSpace( *level, conforming_type ) ) );
        }

        const auto vector_bc =
            std::make_shared<const VectorConformingBasisComplex>( primal.basisComplexPtr(), conforming_type );
        return VectorConformingHierarchicalMultiPatchSplineSpace( vector_bc, vector_levels );
    }

    VectorConformingHierarchicalMultiPatchSplineSpace
        buildHDivHierarchicalMultiPatchSplineSpace( const HierarchicalMultiPatchSplineSpace& primal )
    {
        return buildVectorConformingHierarchicalMultiPatchSplineSpace( primal, ConformingType::Divergence );
    }

    VectorConformingHierarchicalMultiPatchSplineSpace
        buildHCurlHierarchicalMultiPatchSplineSpace( const HierarchicalMultiPatchSplineSpace& primal )
    {
        return buildVectorConformingHierarchicalMultiPatchSplineSpace( primal, ConformingType::Curl );
    }

    HierarchicalMultiPatchSplineSpace
        buildL2HierarchicalMultiPatchSplineSpace( const VectorConformingHierarchicalMultiPatchSplineSpace& hdiv )
    {
        if( hdiv.basisComplex().conformingType() != ConformingType::Divergence )
            throw std::invalid_argument( "Hierarchical L2 construction expects an H(div) hierarchical space." );

        std::vector<std::shared_ptr<const MultiPatchSplineSpace>> l2_levels;
        l2_levels.reserve( hdiv.refinementLevels().size() );
        for( const auto& level : hdiv.refinementLevels() )
        {
            l2_levels.push_back(
                std::make_shared<const MultiPatchSplineSpace>( buildL2MultiPatchSplineSpace( *level ) ) );
        }

        const auto& cmap =
            dynamic_cast<const topology::HierarchicalMultiPatchCombinatorialMap&>(
                hdiv.basisComplex().parametricAtlas().cmap() );
        return buildHierarchicalSplineSpace( l2_levels, topology::leafElements( cmap ) );
    }

    std::vector<VectorHierarchicalLevelFunctionId>
        patchVectorFidsToRefinementLevelFids( const VectorConformingHierarchicalMultiPatchSplineSpace& ss,
                                              const size_t patch_id )
    {
        if( patch_id >= ss.constituents().size() )
            throw std::invalid_argument( "Patch id is outside the hierarchical vector multipatch space." );

        std::vector<VectorHierarchicalLevelFunctionId> out;
        const auto& constituent = *ss.constituents().at( patch_id );
        for( size_t component = 0; component < constituent.scalarBases().size(); component++ )
        {
            const auto& component_basis = *constituent.scalarBases().at( component );
            for( size_t level_ii = 0; level_ii < ss.refinementLevels().size(); level_ii++ )
            {
                const auto& patch_level = *ss.refinementLevels().at( level_ii )->subSpaces().at( patch_id );
                const std::vector<size_t> offsets = componentOffsets( patch_level );
                const auto& patch_level_map = ss.refinementLevels().at( level_ii )->functionIdMap().at( patch_id );

                for( const FunctionId& scalar_fid : component_basis.activeFunctions().at( level_ii ) )
                {
                    const FunctionId patch_vector_fid( scalar_fid.id() + offsets.at( component ) );
                    const auto& [level_global_fid, orientation] = patch_level_map.at( patch_vector_fid.id() );
                    out.push_back( { level_ii, level_global_fid, orientation } );
                }
            }
        }

        if( out.size() != constituent.numFunctions() )
            throw std::runtime_error( "Patch hierarchical vector refinement-level map is incomplete." );
        return out;
    }

    std::vector<std::pair<size_t, FunctionId>>
        vectorFidsToRefinementLevelFids( const VectorConformingHierarchicalMultiPatchSplineSpace& ss )
    {
        std::vector<std::pair<size_t, FunctionId>> out(
            ss.numFunctions(), { std::numeric_limits<size_t>::max(), FunctionId( 0 ) } );
        std::vector<bool> seen( ss.numFunctions(), false );

        for( size_t patch_ii = 0; patch_ii < ss.constituents().size(); patch_ii++ )
        {
            const auto patch_level_fids = patchVectorFidsToRefinementLevelFids( ss, patch_ii );
            const auto& patch_to_global = ss.functionIdMap().at( patch_ii );
            if( patch_level_fids.size() != patch_to_global.size() )
                throw std::runtime_error( "Patch-level and global vector function maps have inconsistent sizes." );

            for( size_t local_fid = 0; local_fid < patch_to_global.size(); local_fid++ )
            {
                const FunctionId global_fid = patch_to_global.at( local_fid ).first;
                const std::pair<size_t, FunctionId> level_fid{
                    patch_level_fids.at( local_fid ).level, patch_level_fids.at( local_fid ).function };

                if( seen.at( global_fid.id() ) )
                {
                    if( out.at( global_fid.id() ).first != level_fid.first or
                        out.at( global_fid.id() ).second.id() != level_fid.second.id() )
                        throw std::runtime_error( "A global hierarchical vector function maps to inconsistent levels." );
                }
                else
                {
                    out.at( global_fid.id() ) = level_fid;
                    seen.at( global_fid.id() ) = true;
                }
            }
        }

        for( const bool was_seen : seen )
        {
            if( not was_seen )
                throw std::runtime_error(
                    "Not all hierarchical vector functions have a corresponding refinement-level function." );
        }
        return out;
    }
}
