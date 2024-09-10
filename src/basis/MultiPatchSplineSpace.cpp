#include <MultiPatchSplineSpace.hpp>

namespace basis
{
    MultiPatchSplineSpace::MultiPatchSplineSpace( const std::shared_ptr<const MultiPatchBasisComplex>& bc,
                                                  const std::vector<std::shared_ptr<const TPSplineSpace>>& constituents,
                                                  const std::vector<std::vector<FunctionId>>& func_ids )
        : mBasisComplex( bc ), mSubSpaces( constituents ), mFuncIds( func_ids )
    {
        if( mFuncIds.size() != mSubSpaces.size() ) throw std::runtime_error( "func_ids and constituents must be the same length" );
        for( size_t i = 0; i < mSubSpaces.size(); i++ )
        {
            if( mFuncIds.at( i ).size() != mSubSpaces.at( i )->numFunctions() )
                throw std::runtime_error( "Bad func_ids input to MultiPatchSplineSpace" );
        }

        mNumFunctions = std::transform_reduce(
            func_ids.begin(), func_ids.end(), Eigen::Index{0}, []( const Eigen::Index accum, const Eigen::Index max_elem ) {
                return std::max( accum, max_elem );
            }, []( const std::vector<FunctionId>& fid ) { return std::max_element( fid.begin(), fid.end() )->id(); } ) + 1;
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
}