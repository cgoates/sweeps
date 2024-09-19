#include <HierarchicalTPSplineSpace.hpp>
#include <set>
#include <algorithm>
#include <Logging.hpp>//FIXME
#include <SparseMatrixUtilities.hpp>
#include <Eigen/SparseLU>

using namespace basis;

std::vector<std::vector<FunctionId>>
    activeFuncs( const topology::HierarchicalTPCombinatorialMap& cmap,
                 const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels )
{
    const std::vector<std::vector<topology::Cell>> leaf_elements = leafElements( cmap );
    std::set<topology::Cell> pruned_cells;
    std::vector<std::vector<FunctionId>> active_funcs;

    for( size_t i = 0; i < refinement_levels.size(); i++ )
    {
        const auto& level_ss = refinement_levels.at( i );
        std::set<topology::Cell> next_pruned_cells;
        std::set<FunctionId> level_active_funcs;
        for( const topology::Cell& leaf_elem : leaf_elements.at( i ) )
        {
            const std::vector<FunctionId> conn = level_ss->connectivity( leaf_elem );
            level_active_funcs.insert( conn.begin(), conn.end() );
            cmap.iterateChildren( leaf_elem, i, [&]( const topology::Cell& c ) {
                next_pruned_cells.insert( c );
                return true;
            } );
        }

        for( const topology::Cell& pruned_cell : pruned_cells )
        {
            const std::vector<FunctionId> conn = level_ss->connectivity( pruned_cell );
            for( const FunctionId& fid : conn ) level_active_funcs.erase( fid );
            cmap.iterateChildren( pruned_cell, i, [&]( const topology::Cell& c ) {
                next_pruned_cells.insert( c );
                return true;
            } );
        }
        pruned_cells = next_pruned_cells;
        active_funcs.emplace_back( level_active_funcs.begin(), level_active_funcs.end() );
    }
    return active_funcs;
}

HierarchicalTPSplineSpace::HierarchicalTPSplineSpace(
    const std::shared_ptr<const HierarchicalTPBasisComplex>& bc,
    const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels )
    : mBasisComplex( bc ),
      mRefinementLevels( refinement_levels ),
      mActiveFunctions( activeFuncs( bc->parametricAtlas().cmap(), refinement_levels ) )
{
    const auto active_mat = [this]( const size_t level_ii ) {
        Eigen::SparseMatrix<double> A( mActiveFunctions.at( level_ii ).size(), mRefinementLevels.at( level_ii )->numFunctions() );
        A.reserve( Eigen::VectorXi::Constant( A.cols(), 1 ) );
        size_t row_ii = 0;
        for( const FunctionId& fid : mActiveFunctions.at( level_ii ) )
            A.coeffRef( row_ii++, fid ) = 1.0;

        return A;
    };
    const auto active_mask = [this]( const size_t level_ii ) {
        const Eigen::Index num_funcs = mRefinementLevels.at( level_ii )->numFunctions();
        Eigen::SparseMatrix<double> A( num_funcs, num_funcs );
        A.reserve( Eigen::VectorXi::Constant( A.cols(), 1 ) );
        for( Eigen::Index col_ii = 0; col_ii < A.cols(); col_ii++ )
            if( std::find( mActiveFunctions.at( level_ii ).begin(), mActiveFunctions.at( level_ii ).end(), FunctionId( col_ii ) ) ==
                mActiveFunctions.at( level_ii ).end() )
                A.coeffRef( col_ii, col_ii ) = 1.0;

        return A;
    };

    mLevelExtractionOps.push_back( active_mat( 0 ) );
    for( size_t i = 1; i < mRefinementLevels.size(); i++ )
    {
        const auto D = refinementOp( *mRefinementLevels.at( i - 1 ), *mRefinementLevels.at( i ), 1e-10 );
        const Eigen::SparseMatrix<double> S = util::verticalConcat( mLevelExtractionOps.back() * D * active_mask( i ), active_mat( i ) );

        mLevelExtractionOps.push_back( S );
    }
}

const HierarchicalTPBasisComplex& HierarchicalTPSplineSpace::basisComplex() const
{
    return *mBasisComplex;
}
const std::shared_ptr<const HierarchicalTPBasisComplex>& HierarchicalTPSplineSpace::basisComplexPtr() const
{
    return mBasisComplex;
}

Eigen::SparseMatrix<double> nonzeroRowsOfColumns( const Eigen::SparseMatrix<double>& mat, const std::vector<FunctionId>& cols )
{
    std::set<Eigen::Index> nonzero_rows;
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve( mat.nonZeros() ); // Overkill, but I don't see a better number to use here.

    for( size_t col_ii = 0; col_ii < cols.size(); ++col_ii )
    {
        for( Eigen::SparseMatrix<double>::InnerIterator it( mat, cols.at( col_ii ) ); it; ++it )
        {
            nonzero_rows.insert( it.row() );
            triplets.emplace_back( it.row(), col_ii, it.value() );
        }
    }

    std::for_each( triplets.begin(),
                   triplets.end(),
                   [&nonzero_rows]( Eigen::Triplet<double>& trip ) {
                       trip = Eigen::Triplet<double>(
                           std::distance( nonzero_rows.begin(), std::ranges::find( nonzero_rows, trip.row() ) ),
                           trip.col(),
                           trip.value() );
                   } );

    Eigen::SparseMatrix<double> result( nonzero_rows.size(), cols.size() );
    result.setFromTriplets( triplets.begin(), triplets.end() );

    return result;
}

Eigen::MatrixXd HierarchicalTPSplineSpace::extractionOperator( const topology::Cell& c ) const
{
    const auto [ level, level_d ] = mBasisComplex->parametricAtlas().cmap().unrefinedAncestorDart( c.dart() );
    const topology::Cell level_c( level_d, c.dim() );
    const std::vector<FunctionId> level_conn = mRefinementLevels.at( level )->connectivity( level_c );
    const Eigen::MatrixXd level_op = mRefinementLevels.at( level )->extractionOperator( level_c );
    return nonzeroRowsOfColumns( mLevelExtractionOps.at( level ), level_conn ) * level_op;
}

std::vector<FunctionId> HierarchicalTPSplineSpace::connectivity( const topology::Cell& c ) const
{
    const auto [ level, level_d ] = mBasisComplex->parametricAtlas().cmap().unrefinedAncestorDart( c.dart() );
    const std::vector<FunctionId> level_conn = mRefinementLevels.at( level )->connectivity( topology::Cell( level_d, c.dim() ) );
    
    std::set<FunctionId> nonzero_rows;
    for( const FunctionId& fid : level_conn )
    {
        for( Eigen::SparseMatrix<double>::InnerIterator it( mLevelExtractionOps.at( level ), fid ); it; ++it )
            nonzero_rows.insert( it.row() );
    }
    return std::vector<FunctionId>( nonzero_rows.begin(), nonzero_rows.end() );
}

size_t HierarchicalTPSplineSpace::numFunctions() const
{
    return mLevelExtractionOps.back().cols();
}

namespace basis
{
    HierarchicalTPSplineSpace
        buildHierarchicalSplineSpace( const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels,
                                      const std::vector<std::vector<topology::Cell>>& leaf_elements )
    {
        std::vector<std::shared_ptr<const TPBasisComplex>> bc_levels;
        bc_levels.reserve( refinement_levels.size() );
        std::vector<std::shared_ptr<const param::TPParametricAtlas>> atlas_levels;
        atlas_levels.reserve( refinement_levels.size() );
        std::vector<std::shared_ptr<const topology::TPCombinatorialMap>> cmap_levels;
        cmap_levels.reserve( refinement_levels.size() );

        std::set<topology::Cell> pruned_cells;

        for( const auto& ss : refinement_levels )
        {
            bc_levels.push_back( ss->basisComplexPtr() );
            atlas_levels.push_back( bc_levels.back()->parametricAtlasPtr() );
            cmap_levels.push_back( atlas_levels.back()->cmapPtr() );
        }

        const auto cmap = std::make_shared<const topology::HierarchicalTPCombinatorialMap>( cmap_levels, leaf_elements );
        const auto atlas = std::make_shared<const param::HierarchicalTPParametricAtlas>( cmap, atlas_levels );
        const auto bc = std::make_shared<const HierarchicalTPBasisComplex>( atlas, bc_levels );

        return HierarchicalTPSplineSpace( bc, refinement_levels ); // NOTE: leaf_elements are recalculated here. Fix if this is a bottleneck.
    }
}