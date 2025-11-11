#include <VectorConformingMultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <MultiPatchSplineSpace.hpp>
#include <GlobalCellMarker.hpp>
#include <ranges>
#include <iostream>
#include <UnionFind.hpp>
#include <Logging.hpp>

using namespace basis;

constexpr bool DEBUG_OUTPUT = false;

const param::MultiPatchParametricAtlas& castAndCheck( const param::ParametricAtlas& atlas )
{
    const auto* mp_atlas = dynamic_cast<const param::MultiPatchParametricAtlas*>( &atlas );
    if( not mp_atlas )
    {
        throw std::runtime_error( "VectorConformingMultiPatchSplineSpace requires a MultiPatchParametricAtlas" );
    }
    return *mp_atlas;
}

bool isCornerCellRelevantForComponent(
    const param::TPParametricAtlas& param,
    const topology::Cell& corner_cell,
    const size_t vector_component,
    const basis::ConformingType conforming_type )
{
    const size_t dim = param.cmap().dim();
    switch( conforming_type )
    {
        case basis::ConformingType::Divergence:
        {
            // For div-conforming, a corner cell is relevant if it lies on the boundary
            // in the direction corresponding to the vector component.
            const param::BaryCoordIsZeroVec corner_bdry = parentDomainBoundary( param, corner_cell );
            return corner_bdry.at( 2 * vector_component ) or corner_bdry.at( 2 * vector_component + 1 );
        }
        case basis::ConformingType::Curl:
        {
            // For curl-conforming, a corner cell is relevant if it lies on any boundary
            // that does not correspond to the vector component.
            const param::BaryCoordIsZeroVec corner_bdry = parentDomainBoundary( param, corner_cell );
            bool relevant = false;
            for( size_t i = 0; i < dim; i++ )
            {
                if( i == vector_component ) continue;
                relevant = relevant or corner_bdry.at( 2 * i ) or corner_bdry.at( 2 * i + 1 );
            }
            return relevant;
        }
        default:
            throw std::runtime_error( "Unknown conforming type" );
    }
}

util::UnionFind unitedFunctions( const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>& constituents,
                                 const param::MultiPatchParametricAtlas& atlas,
                                 const basis::ConformingType& conforming_type,
                                 const std::vector<size_t>& constituent_function_offsets )
{
    const topology::MultiPatchCombinatorialMap& multi_cmap = atlas.cmap();
    const size_t cell_dim = multi_cmap.dim() - 1;

    util::UnionFind uf( constituent_function_offsets.back() );
    for( size_t patch_ii = 0; patch_ii < constituents.size(); patch_ii++ )
    {
        const auto& constituent = *constituents.at( patch_ii );
        const param::TPParametricAtlas& patch_param = *atlas.constituents().at( patch_ii );

        const auto corner_cells = param::cornerCells( patch_param, cell_dim );
        for( const auto corner_cell : corner_cells )
        {
            /// HERE WE COULD THEN do what's inside this for loop for each of the component scalar TP spline spaces
            /// in the constituent but we'll need to check if the given corner cell is supposed to be connected in that
            /// vector component.
            size_t this_component_offset = 0;
            for( size_t vec_comp = 0; vec_comp < constituent.numVectorComponents(); vec_comp++ )
            {
                const auto& constituent_component = *constituent.scalarTPBases().at( vec_comp );
                if( not isCornerCellRelevantForComponent(
                        patch_param, corner_cell, vec_comp, conforming_type ) )
                {
                    LOG( DEBUG_OUTPUT ) << "Skipping patch " << patch_ii << " component " << vec_comp
                                        << " for corner cell " << corner_cell << "\n";
                    this_component_offset += constituent_component.numFunctions();
                    continue;
                }

                const auto [this_lengths, this_order, this_direction] = getIterVars( constituent_component, corner_cell );
                const topology::Cell glob_corner_cell( multi_cmap.toGlobalDart( patch_ii, corner_cell.dart() ), corner_cell.dim() );

                const auto maybe_phi3 = phi( multi_cmap, 3, glob_corner_cell.dart() );
                if( maybe_phi3.has_value() )
                {
                    const topology::Cell glob_neighbor_elem( maybe_phi3.value(), cell_dim );
                    const auto [other_patch_ii, other_corner_cell_dart] = multi_cmap.toLocalDart( glob_neighbor_elem.dart() );
                    
                    constexpr bool reverse = true; // Because phi3 reverses orientation across the face.

                    const auto& other_constituent = *constituents.at( other_patch_ii );
                    /// Here we need to find which vector component corresponds to the vector component from the other patch.
                    const auto [other_vec_comp, aligned] = coordinateTransform( atlas, glob_corner_cell ).at( vec_comp );

                    size_t other_component_offset = 0;
                    for( size_t i = 0; i < other_vec_comp; i++ )
                    {
                        other_component_offset += other_constituent.scalarTPBases().at( i )->numFunctions();
                    }

                    const auto& other_constituent_component = *other_constituent.scalarTPBases().at( other_vec_comp );

                    const auto [other_lengths, other_order, other_direction] = getIterVars(
                        other_constituent_component, topology::Cell( other_corner_cell_dart, cell_dim ), reverse );

                    LOG( DEBUG_OUTPUT ) << "Connecting patch " << patch_ii << " component " << vec_comp
                                        << " with patch " << other_patch_ii << " component " << other_vec_comp << "\n";
                    LOG( DEBUG_OUTPUT ) << "  This lengths: " << this_lengths << ", order: " << this_order
                                        << ", direction: " << this_direction << "\n";
                    LOG( DEBUG_OUTPUT ) << "  Other lengths: " << other_lengths << ", order: " << other_order
                                        << ", direction: " << other_direction << "\n";

                    // Connect the functions from neighboring patches
                    util::iterateTensorProductSynchronized(
                        this_lengths,
                        other_lengths,
                        this_order,
                        other_order,
                        this_direction,
                        other_direction,
                        [&, this_lengths = this_lengths, other_lengths = other_lengths](
                            const util::IndexVec& this_iv, const util::IndexVec& other_iv ) {
                            const FunctionId this_component_fid = util::flatten( this_iv, this_lengths );
                            const FunctionId other_component_fid = util::flatten( other_iv, other_lengths );
                            const FunctionId this_fid = this_component_fid + this_component_offset;
                            const FunctionId other_fid = other_component_fid + other_component_offset;

                            LOG( DEBUG_OUTPUT ) << "    Uniting local func " << this_fid << " with " << other_fid
                                                << " (aligned: " << aligned << ")\n";

                            uf.unite( constituent_function_offsets.at( patch_ii ) + this_fid.id(),
                                        constituent_function_offsets.at( other_patch_ii ) + other_fid.id(), aligned );

                            return true;
                        } );
                }
                this_component_offset += constituent_component.numFunctions();
            }
        }
    }

    return uf;
}

VectorConformingMultiPatchSplineSpace::VectorConformingMultiPatchSplineSpace(
    const std::shared_ptr<const VectorConformingBasisComplex>& bc,
    const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>& constituents )
    : mBasisComplex( bc ), mSubSpaces( constituents ), mParametricAtlas( castAndCheck( bc->parametricAtlas() ) )
{
    std::vector<size_t> constituent_function_offsets;
    constituent_function_offsets.push_back( 0 );
    for( size_t i = 0; i < constituents.size(); i++ )
    {
        constituent_function_offsets.push_back( constituent_function_offsets.back() + constituents.at( i )->numFunctions() );
    }

    util::UnionFind uf = unitedFunctions( constituents, mParametricAtlas, bc->conformingType(), constituent_function_offsets );

    mNumFunctions = uf.numSets();
    mFuncIds.reserve( constituents.size() );

    // Pull out the function ids from the union find
    std::unordered_map<size_t, FunctionId> fid_map;
    for( size_t patch_ii = 0; patch_ii < constituents.size(); patch_ii++ )
    {
        const size_t offset = constituent_function_offsets.at( patch_ii );
        mFuncIds.push_back( std::vector<std::pair<FunctionId, bool>>() );
        std::vector<std::pair<FunctionId, bool>>& constituent_fids = mFuncIds.back();
        constituent_fids.reserve( constituents.at( patch_ii )->numFunctions() );
        for( size_t func_ii = 0, num_funcs = constituents.at( patch_ii )->numFunctions(); func_ii < num_funcs; func_ii++ )
        {
            const auto pr = uf.findWithOrientation( offset + func_ii );
            const auto it = fid_map.find( pr.first );
            if( it == fid_map.end() )
            {
                const FunctionId new_fid( fid_map.size() );
                fid_map.emplace( pr.first, new_fid );
                constituent_fids.push_back( { new_fid, pr.second } );
            }
            else
            {
                constituent_fids.push_back( { fid_map.at( pr.first ), pr.second } );
            }
        }
    }
}

const VectorConformingBasisComplex& VectorConformingMultiPatchSplineSpace::basisComplex() const
{
    return *mBasisComplex;
}

size_t VectorConformingMultiPatchSplineSpace::numVectorComponents() const
{
    return mParametricAtlas.cmap().dim();
}

Eigen::MatrixXd VectorConformingMultiPatchSplineSpace::extractionOperator( const topology::Cell& c ) const
{
    const auto [patch_id, local_d] = mParametricAtlas.cmap().toLocalDart( c.dart() );
    const topology::Cell local_c( local_d, c.dim() );
    const std::vector<FunctionId> local_conn = mSubSpaces.at( patch_id )->connectivity( local_c );

    const auto patch_func_map = mFuncIds.at( patch_id );

    Eigen::VectorXd orientations( local_conn.size() );
    for( size_t i = 0; i < local_conn.size(); i++ )
    {
        const auto& [global_fid, orientation] = patch_func_map.at( local_conn.at( i ).id() );
        orientations(i) = orientation ? 1.0 : -1.0;
    }
    return orientations.asDiagonal() * mSubSpaces.at( patch_id )->extractionOperator( local_c );
}

std::vector<FunctionId> VectorConformingMultiPatchSplineSpace::connectivity( const topology::Cell& c ) const
{
    const auto [patch_id, local_d] = mParametricAtlas.cmap().toLocalDart( c.dart() );
    const topology::Cell local_c( local_d, c.dim() );
    const std::vector<FunctionId> local_conn = mSubSpaces.at( patch_id )->connectivity( local_c );

    const auto& patch_func_map = mFuncIds.at( patch_id );

    std::vector<FunctionId> global_conn;
    global_conn.reserve( local_conn.size() );
    std::transform( local_conn.begin(),
                    local_conn.end(),
                    std::back_inserter( global_conn ),
                    [&]( const FunctionId& local_fid ) { return patch_func_map.at( local_fid.id() ).first; } );

    return global_conn;
}

size_t VectorConformingMultiPatchSplineSpace::numFunctions() const
{
    return mNumFunctions;
}