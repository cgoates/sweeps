#include <MultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>
#include <ranges>
#include <GlobalCellMarker.hpp>

namespace basis
{
    std::tuple<util::IndexVec, util::IndexVec, SmallVector<std::variant<bool, size_t>, 3>>
        getIterVars( const TPSplineSpace& constituent, const topology::Cell& corner, const bool reverse_dart = false )
    {
        const auto& param = constituent.basisComplex().parametricAtlas();
        const size_t param_dim = param.cmap().dim();
        const util::IndexVec lengths = getTPLengths( constituent );
        const param::ParentDomain pd = param.parentDomain( topology::Cell( corner.dart(), param_dim ) );
        const param::BaryCoordIsZeroVec corner_bdry = param::parentDomainBoundary( param, corner );
        SmallVector<std::variant<bool, size_t>, 3> direction;
        util::IndexVec order( param_dim );
        switch( corner.dim() )
        {
            case 0:
            {
                const param::BaryCoordIsZeroVec corner_bdry = param::parentDomainBoundary( constituent.basisComplex().parametricAtlas(), corner );
                iterateGroups( pd, [&]( const size_t first_idx, const auto, const auto ){
                    if( corner_bdry.at( first_idx ) )
                        direction.push_back( lengths.at( direction.size() ) - 1 );
                    else if( corner_bdry.at( first_idx + 1 ) )
                        direction.push_back( size_t{0} );
                    else
                        throw std::runtime_error( "Bad boundary for a vertex" );
                } );
                return { lengths, std::ranges::to<util::IndexVec>( std::ranges::iota_view( size_t{0}, param_dim ) ), direction };
            }
            case 1:
            {
                size_t dummy_index = 1;
                const param::BaryCoordIsZeroVec corner_bdry = param::parentDomainBoundary( constituent.basisComplex().parametricAtlas(), corner );
                iterateGroups( pd, [&]( const size_t first_idx, const auto, const auto ){
                    if( corner_bdry.at( first_idx ) )
                    {
                        order.at( dummy_index++ ) = direction.size();
                        direction.push_back( lengths.at( direction.size() ) - 1 );
                    }
                    else if( corner_bdry.at( first_idx + 1 ) )
                    {
                        order.at( dummy_index++ ) = direction.size();
                        direction.push_back( size_t{0} );
                    }
                    else
                    {
                        order.at( 0 ) = direction.size();
                        const param::ParentPoint ppt = param.parentPoint( corner.dart() );
                        // If the vertex has compressed coordinate zero, then we're going in a positive direction unless it's reversed.
                        direction.push_back( reverse_dart != ppt.mBaryCoordIsZero.at( first_idx + 1 ) );
                    }
                } );
                return { lengths, order, direction };
            }
            break;
            case 2:
            {
                // Build a coordinate frame
                const auto [ ppt00, ppt10, ppt01 ] = [&]() -> std::tuple<param::ParentPoint, param::ParentPoint, param::ParentPoint> {
                    const param::ParentPoint ppt00 = param.parentPoint( corner.dart() );
                    const param::ParentPoint ppt10 = param.parentPoint( phi( param.cmap(), 1, corner.dart() ).value() );
                    if( reverse_dart )
                    {
                        const param::ParentPoint ppt01 =
                            param.parentPoint( phi( param.cmap(), { 1, 1 }, corner.dart() ).value() );
                        return { ppt10, ppt00, ppt01 };
                    }
                    else
                    {
                        const param::ParentPoint ppt01 =
                            param.parentPoint( phi( param.cmap(), -1, corner.dart() ).value() );
                        return { ppt00, ppt10, ppt01 };
                    }
                }();

                const param::BaryCoordIsZeroVec corner_bdry = param::join( param::join( ppt00.mBaryCoordIsZero, ppt10.mBaryCoordIsZero ), ppt01.mBaryCoordIsZero );
                iterateGroups( pd, [&, ppt00 = ppt00, ppt10 = ppt10, ppt01 = ppt01]( const size_t first_idx, const auto, const auto ){
                    if( corner_bdry.at( first_idx ) )
                    {
                        order.at( 2 ) = direction.size();
                        direction.push_back( lengths.at( direction.size() ) - 1 );
                    }
                    else if( corner_bdry.at( first_idx + 1 ) )
                    {
                        order.at( 2 ) = direction.size();
                        direction.push_back( size_t{0} );
                    }
                    else if( ppt10.mBaryCoordIsZero.at( first_idx ) and not ppt00.mBaryCoordIsZero.at( first_idx ) )
                    {
                        order.at( 0 ) = direction.size();
                        direction.push_back( true );
                    }
                    else if( not ppt10.mBaryCoordIsZero.at( first_idx ) and ppt00.mBaryCoordIsZero.at( first_idx ) )
                    {
                        order.at( 0 ) = direction.size();
                        direction.push_back( false );
                    }
                    else if( ppt01.mBaryCoordIsZero.at( first_idx ) and not ppt00.mBaryCoordIsZero.at( first_idx ) )
                    {
                        order.at( 1 ) = direction.size();
                        direction.push_back( true );
                    }
                    else if( not ppt01.mBaryCoordIsZero.at( first_idx ) and ppt00.mBaryCoordIsZero.at( first_idx ) )
                    {
                        order.at( 1 ) = direction.size();
                        direction.push_back( false );
                    }
                } );

                return { lengths, order, direction };
            }
            break;
            default:
                throw std::runtime_error( "Corner cells should be dimension 2 or less" );
        }
    }

    bool onSameVertex( const topology::MultiPatchCombinatorialMap& cmap, const topology::Dart& d1, const topology::Dart& d2 )
    {
        topology::LocalCellMarker m( 0 );
        m.mark( cmap, topology::Vertex( d1 ) );
        return m.isMarked( topology::Vertex( d2 ) );
    }

    MultiPatchSplineSpace::MultiPatchSplineSpace( const std::shared_ptr<const MultiPatchBasisComplex>& bc,
                                                  const std::vector<std::shared_ptr<const TPSplineSpace>>& constituents )
        : mBasisComplex( bc ), mSubSpaces( constituents )
    {
        const topology::MultiPatchCombinatorialMap& multi_cmap = bc->parametricAtlas().cmap();
        const size_t param_dim = multi_cmap.dim();
        mFuncIds.reserve( constituents.size() );
        for( size_t i = 0; i < constituents.size(); i++ )
        {
            mFuncIds.push_back( std::vector<FunctionId>( constituents.at( i )->numFunctions(), -1 ) );
        }

        mNumFunctions = 0;
        for( size_t patch_ii = 0; patch_ii < constituents.size(); patch_ii++ )
        {
            const TPSplineSpace& constituent = *constituents.at( patch_ii );
            auto& constituent_fids = mFuncIds.at( patch_ii );
            // Fill every fid that hasn't been set already with the next available ids
            for( size_t func_ii = 0; func_ii < constituent.numFunctions(); func_ii++ )
            {
                if( constituent_fids.at( func_ii ) == -1 )
                {
                    constituent_fids.at( func_ii ) = mNumFunctions;
                    mNumFunctions = mNumFunctions + 1;
                }
            }

            // Iterate the corner cells to connect to neighboring patches
            for( size_t cell_dim = 0; cell_dim < param_dim; cell_dim++ )
            {
                const auto corner_cells = param::cornerCells( constituent.basisComplex().parametricAtlas(), cell_dim );
                for( const auto corner_cell : corner_cells )
                {
                    const auto [this_lengths, this_order, this_direction] = getIterVars( constituent, corner_cell );
                    const topology::Cell glob_corner_cell( multi_cmap.toGlobalDart( patch_ii, corner_cell.dart() ), corner_cell.dim() );

                    // Iterate adjacent patches
                    iterateAdjacentCells( multi_cmap, glob_corner_cell, param_dim, [&]( const topology::Cell& glob_neighbor_elem ) {
                        const auto [other_patch_ii, other_corner_cell_dart] = multi_cmap.toLocalDart( glob_neighbor_elem.dart() );
                        if( other_patch_ii == patch_ii ) return true; // No need to connect to itself

                        const bool reverse =
                            not onSameVertex( multi_cmap, glob_corner_cell.dart(), glob_neighbor_elem.dart() );

                        const auto [other_lengths, other_order, other_direction] =
                            getIterVars( *constituents.at( other_patch_ii ),
                                         topology::Cell( other_corner_cell_dart, cell_dim ),
                                         reverse );

                        std::vector<FunctionId>& other_constituent_fids = mFuncIds.at( other_patch_ii );

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
                                const FunctionId this_fid = util::flatten( this_iv, this_lengths );
                                const FunctionId other_fid = util::flatten( other_iv, other_lengths );

                                if( constituent_fids.at( this_fid ) == -1 or
                                    ( other_constituent_fids.at( other_fid ) != -1 and
                                      other_constituent_fids.at( other_fid ) != constituent_fids.at( this_fid ) ) )
                                {
                                    std::cerr << constituent_fids.at( this_fid ) << " vs " << other_constituent_fids.at( other_fid ) << std::endl;
                                    throw std::runtime_error( "Problem connecting neighboring functions!" );
                                }

                                other_constituent_fids.at( other_fid ) = constituent_fids.at( this_fid );

                                return true;
                            } );

                        return true;
                    } );
                }
            }
        }
    }

    const MultiPatchBasisComplex& MultiPatchSplineSpace::basisComplex() const
    {
        return *mBasisComplex;
    }

    Eigen::MatrixXd MultiPatchSplineSpace::extractionOperator( const topology::Cell& c ) const
    {
        const auto [patch_id, local_d] = mBasisComplex->parametricAtlas().cmap().toLocalDart( c.dart() );
        const topology::Cell local_c( local_d, c.dim() );
        return mSubSpaces.at( patch_id )->extractionOperator( local_c );
    }

    std::vector<FunctionId> MultiPatchSplineSpace::connectivity( const topology::Cell& c ) const
    {
        const auto [patch_id, local_d] = mBasisComplex->parametricAtlas().cmap().toLocalDart( c.dart() );
        const topology::Cell local_c( local_d, c.dim() );
        const std::vector<FunctionId> local_conn = mSubSpaces.at( patch_id )->connectivity( local_c );

        const auto& patch_func_map = mFuncIds.at( patch_id );

        std::vector<FunctionId> global_conn;
        global_conn.reserve( local_conn.size() );
        std::transform( local_conn.begin(),
                        local_conn.end(),
                        std::back_inserter( global_conn ),
                        [&]( const FunctionId& local_fid ) { return patch_func_map.at( local_fid.id() ); } );

        return global_conn;
    }

    size_t MultiPatchSplineSpace::numFunctions() const
    {
        return mNumFunctions;
    }

    MultiPatchSplineSpace buildMultiPatchSplineSpace(
        const std::vector<std::shared_ptr<const TPSplineSpace>>& patches,
        const std::map<std::pair<size_t, topology::Dart>, std::pair<size_t, topology::Dart>>& connections )
    {
        std::vector<std::shared_ptr<const TPBasisComplex>> bc_patches;
        bc_patches.reserve( patches.size() );
        std::vector<std::shared_ptr<const param::TPParametricAtlas>> atlas_patches;
        atlas_patches.reserve( patches.size() );
        std::vector<std::shared_ptr<const topology::TPCombinatorialMap>> cmap_patches;
        cmap_patches.reserve( patches.size() );

        std::set<topology::Cell> pruned_cells;

        for( const auto& ss : patches )
        {
            bc_patches.push_back( ss->basisComplexPtr() );
            atlas_patches.push_back( bc_patches.back()->parametricAtlasPtr() );
            cmap_patches.push_back( atlas_patches.back()->cmapPtr() );
        }

        const auto cmap = std::make_shared<const topology::MultiPatchCombinatorialMap>( cmap_patches, connections );
        const auto atlas = std::make_shared<const param::MultiPatchParametricAtlas>( cmap, atlas_patches );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>( atlas, bc_patches );

        return MultiPatchSplineSpace( bc, patches );
    }
}