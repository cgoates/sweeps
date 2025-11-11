#include <MultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>
#include <ranges>
#include <GlobalCellMarker.hpp>
#include <iostream>

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

        const auto handle_on_group_boundary = [&]( const size_t first_idx, const size_t order_idx ) {
            if( corner_bdry.at( first_idx ) )
            {
                order.at( order_idx ) = direction.size();
                direction.push_back( lengths.at( direction.size() ) - 1 );
                return true;
            }
            else if( corner_bdry.at( first_idx + 1 ) )
            {
                order.at( order_idx ) = direction.size();
                direction.push_back( size_t{0} );
                return true;
            }

            return false;
        };

        switch( corner.dim() )
        {
            case 0:
            {
                iterateGroups( pd, [&]( const size_t first_idx, const auto, const auto ){
                    if( not handle_on_group_boundary( first_idx, 0 ) )
                        throw std::runtime_error( "Bad boundary for a vertex" );
                } );
                const auto range = std::ranges::iota_view( size_t{0}, param_dim );
                return { lengths, util::IndexVec( range.begin(), range.end() ), direction };
            }
            case 1:
            {
                size_t next_non_edge_group_id = 1;
                iterateGroups( pd, [&]( const size_t first_idx, const auto, const auto ){
                    if( handle_on_group_boundary( first_idx, next_non_edge_group_id ) )
                        next_non_edge_group_id++;
                    else // This is the group along which the edge runs
                    {
                        order.at( 0 ) = direction.size();
                        const param::ParentPoint ppt = param.parentPoint( corner.dart() );
                        // If the vertex has compressed coordinate zero, then we're going in a positive direction unless it's reversed.
                        direction.push_back( reverse_dart != ppt.mBaryCoordIsZero.at( first_idx + 1 ) );
                    }
                } );
                return { lengths, order, direction };
            }
            case 2:
            {
                // Build a coordinate frame
                const SmallVector<param::ParentPoint, 4> frame_vec = param::getFrame( param, corner, reverse_dart );

                const auto find_face_direction = [&]( const param::ParentPoint& ppt_a,
                                                      const param::ParentPoint& ppt_b,
                                                      const size_t first_idx,
                                                      const size_t order_idx ) {
                    if( ppt_b.mBaryCoordIsZero.at( first_idx ) and not ppt_a.mBaryCoordIsZero.at( first_idx ) )
                    {
                        order.at( order_idx ) = direction.size();
                        direction.push_back( true );
                        return true;
                    }
                    else if( not ppt_b.mBaryCoordIsZero.at( first_idx ) and ppt_a.mBaryCoordIsZero.at( first_idx ) )
                    {
                        order.at( order_idx ) = direction.size();
                        direction.push_back( false );
                        return true;
                    }
                    return false;
                };

                iterateGroups( pd, [&]( const size_t first_idx, const auto, const auto ){
                    if( not handle_on_group_boundary( first_idx, 2 ) )
                    {
                        // Then the face runs along this group.
                        // Figure out which of the two face directions this group corresponds to.
                        if( not find_face_direction( frame_vec.at( 0 ), frame_vec.at( 1 ), first_idx, 0 ) )
                            find_face_direction( frame_vec.at( 0 ), frame_vec.at( 2 ), first_idx, 1 );
                    }
                } );

                return { lengths, order, direction };
            }
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

                                if( constituent_fids.at( this_fid ) == -1 )
                                {
                                    constituent_fids.at( this_fid ) = mNumFunctions;
                                    mNumFunctions = mNumFunctions + 1;
                                }

                                if( other_constituent_fids.at( other_fid ) != -1 and
                                    other_constituent_fids.at( other_fid ) != constituent_fids.at( this_fid ) )
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

            // Fill every fid that hasn't been set already with the next available ids
            for( size_t func_ii = 0, num_funcs = constituent.numFunctions(); func_ii < num_funcs; func_ii++ )
            {
                if( constituent_fids.at( func_ii ) == -1 )
                {
                    constituent_fids.at( func_ii ) = mNumFunctions;
                    mNumFunctions = mNumFunctions + 1;
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

    MultiPatchSplineSpace buildMultiPatchSplineSpace(
        const std::vector<std::shared_ptr<const TPSplineSpace>>& patches,
        const topology::MultiPatchCombinatorialMap::InternalConnectionsMap& connections )
    {
        std::vector<std::shared_ptr<const TPBasisComplex>> bc_patches;
        bc_patches.reserve( patches.size() );
        std::vector<std::shared_ptr<const param::TPParametricAtlas>> atlas_patches;
        atlas_patches.reserve( patches.size() );
        std::vector<std::shared_ptr<const topology::TPCombinatorialMap>> cmap_patches;
        cmap_patches.reserve( patches.size() );

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

    MultiPatchSplineSpace degreeRefineOrCoarsen( const MultiPatchSplineSpace& ss,
                                                 const std::function<DegreeAndKnotVector( const size_t )>& degree_and_kv_func )
    {
        std::vector<std::shared_ptr<const basis::TPSplineSpace>> constituents;
        constituents.reserve( ss.subSpaces().size() );

        for( size_t patch_ii = 0; patch_ii < ss.subSpaces().size(); patch_ii++ )
        {
            const auto new_degrees_and_kvs = degree_and_kv_func( patch_ii );
            constituents.push_back( std::make_shared<const basis::TPSplineSpace>( basis::buildBSpline( new_degrees_and_kvs.kvs, new_degrees_and_kvs.degrees ) ) );
        }
        const auto temp = buildMultiPatchSplineSpace( constituents, ss.basisComplex().parametricAtlas().cmap().connections() );

        const auto bc = std::make_shared<const basis::MultiPatchBasisComplex>( ss.basisComplex().parametricAtlasPtr(), temp.basisComplex().constituents() );
        return basis::MultiPatchSplineSpace( bc, temp.subSpaces() );
    }

    Eigen::MatrixXd multiPatchCoefficients( const MultiPatchSplineSpace& ss,
                                            const std::vector<Eigen::MatrixXd>& patch_coeffs )
    {
        const auto& func_ids = ss.functionIdMap();
        if( patch_coeffs.size() != func_ids.size() )
            throw std::invalid_argument( "Wrong number of patch coefficients in multiPatchCoefficients" );

        Eigen::MatrixXd out = Eigen::MatrixXd::Zero( ss.numFunctions(), patch_coeffs.at( 0 ).cols() );
        std::vector<size_t> num_funcs( ss.numFunctions(), 0 );

        for( size_t patch_ii = 0; patch_ii < patch_coeffs.size(); patch_ii++ )
        {
            const auto& patch_coeff = patch_coeffs.at( patch_ii );
            const auto& patch_func_ids = func_ids.at( patch_ii );
            if( patch_coeff.rows() != (Eigen::Index)patch_func_ids.size() )
                throw std::invalid_argument( "Wrong number of control points for patch " + std::to_string( patch_ii ) + " in multiPatchCoefficients. "
                                              "Expected " + std::to_string( patch_func_ids.size() ) + ", got " + std::to_string( patch_coeff.rows() ) + "." );

            for( size_t func_ii = 0; func_ii < patch_func_ids.size(); func_ii++ )
            {
                const FunctionId global_fid = patch_func_ids.at( func_ii );
                out.row( global_fid ) += patch_coeff.row( func_ii );
                num_funcs.at( global_fid )++;
            }
        }

        for( size_t global_fid = 0; global_fid < num_funcs.size(); global_fid++ )
        {
            if( num_funcs.at( global_fid ) > 1 )
                out.row( global_fid ) /= num_funcs.at( global_fid );
        }

        return out;
    }
}