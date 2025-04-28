#include <HierarchicalMultiPatchSplineSpace.hpp>
#include <set>
#include <CombinatorialMapMethods.hpp>

using namespace basis;

struct ActiveFuncs
{
    std::vector<std::vector<std::vector<FunctionId>>> active_funcs;
    std::vector<std::vector<FunctionId>> patch_level_active_funcs_to_mp_funcs;
    size_t num_active_funcs;
};

ActiveFuncs activeFuncs( const topology::HierarchicalMultiPatchCombinatorialMap& cmap,
                         const std::vector<std::shared_ptr<const MultiPatchSplineSpace>>& refinement_levels )
{
    const std::vector<std::vector<topology::Cell>> leaf_elements = leafElements( cmap );
    std::set<topology::Cell> pruned_cells;
    std::vector<std::vector<std::vector<FunctionId>>> active_funcs(
        refinement_levels.front()->subSpaces().size(),
        std::vector<std::vector<FunctionId>>( refinement_levels.size() ) );

    std::vector<std::vector<FunctionId>> patch_active_to_mp_active( refinement_levels.front()->subSpaces().size() );

    size_t n_active_funcs = 0;

    for( size_t level_ii = 0; level_ii < refinement_levels.size(); level_ii++ )
    {
        const auto& level_ss = refinement_levels.at( level_ii );
        const auto& level_cmap = *cmap.refinementLevels().at( level_ii );
        std::set<topology::Cell> next_pruned_cells;
        std::set<FunctionId> level_active_funcs;
        for( const topology::Cell& leaf_elem : leaf_elements.at( level_ii ) )
        {
            const std::vector<FunctionId> conn = level_ss->connectivity( leaf_elem );
            level_active_funcs.insert( conn.begin(), conn.end() );
            iterateAdjacentCells( level_cmap, leaf_elem, cmap.dim() - 1, [&]( const topology::Cell& c ) {
                const auto maybe_phi = phi( level_cmap, cmap.dim(), c.dart() );
                if( maybe_phi )
                {
                    bool ancestor_leaf = false;
                    if( not ancestor_leaf )
                    {
                        cmap.iterateAncestors( cmap.toGlobalDart( level_ii, maybe_phi.value() ), [&]( const topology::Dart& ancestor ) {
                            if( cmap.isUnrefinedLeafDart( ancestor ) )
                            {
                                ancestor_leaf = true;
                            }
                            return not ancestor_leaf;
                        } );
                    }
                    if( ancestor_leaf )
                    {
                        const std::vector<FunctionId> conn = level_ss->connectivity( topology::Cell( maybe_phi.value(), cmap.dim() ) );
                        for( const FunctionId& fid : conn ) level_active_funcs.erase( fid );
                    }
                }
                return true;
            } );
        }

        for( size_t patch_ii = 0; patch_ii < level_ss->subSpaces().size(); patch_ii++ )
        {
            auto& patch_level_active_funcs = active_funcs.at( patch_ii ).at( level_ii );
            const auto& patch_funcs_to_mp_funcs = level_ss->functionIdMap().at( patch_ii );
            for( size_t patch_fid = 0; patch_fid < patch_funcs_to_mp_funcs.size(); patch_fid++ )
            {
                const auto it = level_active_funcs.find( patch_funcs_to_mp_funcs.at( patch_fid ) );
                if( it != level_active_funcs.end() )
                {
                    patch_level_active_funcs.push_back( patch_fid );
                    patch_active_to_mp_active.at( patch_ii ).push_back( n_active_funcs + std::distance( level_active_funcs.begin(), it ) );
                }
            }
        }
        n_active_funcs += level_active_funcs.size();
    }
    return { active_funcs, patch_active_to_mp_active, n_active_funcs };
}

std::vector<std::shared_ptr<const HierarchicalTPSplineSpace>>
    initializeConstituents( const HierarchicalMultiPatchBasisComplex& bc,
                            const std::vector<std::shared_ptr<const MultiPatchSplineSpace>>& refinement_levels,
                            const std::vector<std::vector<std::vector<FunctionId>>>& active_funcs )
{
    if( bc.refinementLevels().size() != refinement_levels.size() )
        throw std::invalid_argument( "Inconsistent refinement levels and cmap provided to HierarchicalMultiPatchSplineSpace." );

    std::vector<std::vector<std::shared_ptr<const TPSplineSpace>>> constituent_levels( bc.constituents().size() );
    for( const auto& level : refinement_levels )
    {
        if( level->subSpaces().size() != bc.constituents().size() )
            throw std::invalid_argument( "All refinement levels must have the same number of constituents." );
        for( size_t patch_ii = 0; patch_ii < bc.constituents().size(); patch_ii++ )
        {
            constituent_levels.at( patch_ii ).push_back( level->subSpaces().at( patch_ii ) );
        }
    }

    std::vector<std::shared_ptr<const HierarchicalTPSplineSpace>> out;
    out.reserve( bc.constituents().size() );
    for( size_t i = 0; i < bc.constituents().size(); i++ )
    {
        out.push_back( std::make_shared<HierarchicalTPSplineSpace>( bc.constituents().at( i ), constituent_levels.at( i ), active_funcs.at( i ) ) );
    }
    return out;
}

HierarchicalMultiPatchSplineSpace::HierarchicalMultiPatchSplineSpace(
    const std::shared_ptr<const HierarchicalMultiPatchBasisComplex>& bc,
    const std::vector<std::shared_ptr<const MultiPatchSplineSpace>>& refinement_levels )
    : mBasisComplex( bc ), mRefinementLevels( refinement_levels )
{
    const auto active_funcs_struct = activeFuncs( bc->parametricAtlas().cmap(), refinement_levels );
    mConstituents = initializeConstituents( *bc, refinement_levels, active_funcs_struct.active_funcs );
    mNumActiveFuncs = active_funcs_struct.num_active_funcs;
    mFuncIds = active_funcs_struct.patch_level_active_funcs_to_mp_funcs;
}

const HierarchicalMultiPatchBasisComplex& HierarchicalMultiPatchSplineSpace::basisComplex() const { return *mBasisComplex; }
const std::shared_ptr<const HierarchicalMultiPatchBasisComplex>& HierarchicalMultiPatchSplineSpace::basisComplexPtr() const
{
    return mBasisComplex;
}

Eigen::MatrixXd HierarchicalMultiPatchSplineSpace::extractionOperator( const topology::Cell& c ) const
{
    const auto [patch_ii, constituent_d] = mBasisComplex->parametricAtlas().cmap().dartRanges().toLocalDart( c.dart() );
    const topology::Cell constituent_c( constituent_d, c.dim() );
    return mConstituents.at( patch_ii )->extractionOperator( constituent_c );
}

std::vector<FunctionId> HierarchicalMultiPatchSplineSpace::connectivity( const topology::Cell& c ) const
{
    const auto [patch_ii, constituent_d] = mBasisComplex->parametricAtlas().cmap().dartRanges().toLocalDart( c.dart() );
    const topology::Cell constituent_c( constituent_d, c.dim() );
    const std::vector<FunctionId> constituent_conn = mConstituents.at( patch_ii )->connectivity( constituent_c );

    std::vector<FunctionId> conn( constituent_conn.size(), 0 );
    std::transform( constituent_conn.begin(), constituent_conn.end(), conn.begin(), [&]( const FunctionId& fid ) {
        return mFuncIds.at( patch_ii ).at( fid );
    } );
    return conn;
}

size_t HierarchicalMultiPatchSplineSpace::numFunctions() const
{
    return mNumActiveFuncs;
}

namespace basis
{
    HierarchicalMultiPatchSplineSpace
        buildHierarchicalSplineSpace( const std::vector<std::shared_ptr<const MultiPatchSplineSpace>>& refinement_levels,
                                      const std::vector<std::vector<topology::Cell>>& leaf_elements )
    {
        std::vector<std::shared_ptr<const MultiPatchBasisComplex>> bc_levels;
        bc_levels.reserve( refinement_levels.size() );
        std::vector<std::shared_ptr<const param::MultiPatchParametricAtlas>> atlas_levels;
        atlas_levels.reserve( refinement_levels.size() );
        std::vector<std::shared_ptr<const topology::MultiPatchCombinatorialMap>> cmap_levels;
        cmap_levels.reserve( refinement_levels.size() );

        std::set<topology::Cell> pruned_cells;

        for( const auto& ss : refinement_levels )
        {
            bc_levels.push_back( ss->basisComplexPtr() );
            atlas_levels.push_back( bc_levels.back()->parametricAtlasPtr() );
            cmap_levels.push_back( atlas_levels.back()->cmapPtr() );
        }

        const auto cmap = std::make_shared<const topology::HierarchicalMultiPatchCombinatorialMap>( cmap_levels, leaf_elements );
        const auto atlas = std::make_shared<const param::HierarchicalMultiPatchParametricAtlas>( cmap, atlas_levels );
        const auto bc = std::make_shared<const HierarchicalMultiPatchBasisComplex>( atlas, bc_levels );

        return HierarchicalMultiPatchSplineSpace( bc, refinement_levels ); // NOTE: leaf_elements are recalculated here. Fix if this is a bottleneck.
    }
}