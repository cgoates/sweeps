#include <HierarchicalTPSplineSpace.hpp>
#include <set>
#include <algorithm>
#include <SparseMatrixUtilities.hpp>
#include <Eigen/SparseLU>
#include <CombinatorialMapMethods.hpp>

using namespace basis;

std::vector<std::vector<FunctionId>>
    activeFuncs( const topology::HierarchicalTPCombinatorialMap& cmap,
                 const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels )
{
    const std::vector<std::vector<topology::Cell>> leaf_elements = leafElements( cmap );
    std::vector<std::vector<FunctionId>> active_funcs;

    for( size_t i = 0; i < refinement_levels.size(); i++ )
    {
        const auto& level_ss = refinement_levels.at( i );
        const auto& level_cmap = *cmap.refinementLevels().at( i );
        std::set<topology::Cell> cells_to_prune;
        std::set<FunctionId> level_active_funcs;
        for( const topology::Cell& leaf_elem : leaf_elements.at( i ) )
        {
            const std::vector<FunctionId> conn = level_ss->connectivity( leaf_elem );
            level_active_funcs.insert( conn.begin(), conn.end() );
            iterateAdjacentCells( level_cmap, leaf_elem, cmap.dim() - 1, [&]( const topology::Cell& c ) {
                const auto maybe_phi = phi( level_cmap, cmap.dim(), c.dart() );
                if( maybe_phi )
                {
                    bool ancestor_leaf = false;
                    iterateDartsOfCell( level_cmap, topology::Cell( maybe_phi.value(), cmap.dim() ), [&]( const topology::Dart& d ) {
                        cmap.iterateAncestors( cmap.dartRanges().toGlobalDart( i, d ), [&]( const topology::Dart& ancestor ) {
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
        for( const topology::Cell& c : cells_to_prune )
        {
            const std::vector<FunctionId> conn = level_ss->connectivity( c );
            for( const FunctionId& fid : conn ) level_active_funcs.erase( fid );
        }

        active_funcs.emplace_back( level_active_funcs.begin(), level_active_funcs.end() );
    }
    return active_funcs;
}

HierarchicalTPSplineSpace::HierarchicalTPSplineSpace(
    const std::shared_ptr<const HierarchicalTPBasisComplex>& bc,
    const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels )
    : HierarchicalTPSplineSpace( bc, refinement_levels, activeFuncs( bc->parametricAtlas().cmap(), refinement_levels ) )
{}

HierarchicalTPSplineSpace::HierarchicalTPSplineSpace(
    const std::shared_ptr<const HierarchicalTPBasisComplex>& bc,
    const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels,
    const std::vector<std::vector<FunctionId>>& active_funcs )
    : mBasisComplex( bc ),
      mRefinementLevels( refinement_levels ),
      mActiveFunctions( active_funcs )
{
    const auto active_mat = [this]( const size_t level_ii ) {
        Eigen::SparseMatrix<double> A( mActiveFunctions.at( level_ii ).size(), mRefinementLevels.at( level_ii )->numFunctions() );
        A.reserve( Eigen::VectorXi::Constant( A.cols(), 1 ) );
        size_t row_ii = 0;
        for( const FunctionId& fid : mActiveFunctions.at( level_ii ) )
            A.coeffRef( row_ii++, fid ) = 1.0;

        return A;
    };
    const auto active_mask = [this]( const size_t level_ii, Eigen::SparseMatrix<double>& M ) {
        for( const Eigen::Index col_ii : mActiveFunctions.at( level_ii ) )
            M.col( col_ii ) = Eigen::SparseMatrix<double>( M.rows(), 1 );
        M.makeCompressed();
    };

    mLevelExtractionOps.push_back( active_mat( 0 ) );
    for( size_t i = 1; i < mRefinementLevels.size(); i++ )
    {
        // TODO: this can be optimized by only assembling the columns of the refinementOp that are for active functions,
        // since the kronecker product is a slow point for bigger spline spaces.
        // There are also more performant ways to do a product of a kronecker product with another matrix.  Look into that.
        Eigen::SparseMatrix<double> D = refinementOp( *mRefinementLevels.at( i - 1 ), *mRefinementLevels.at( i ), 1e-10 );
        active_mask( i, D );
        Eigen::SparseMatrix<double> temp( mLevelExtractionOps.back().rows(), D.cols() );
        // This is temp = mLevelExtractionOps.back() * D, but using a more performant version, since this is a bottleneck.
        Eigen::internal::conservative_sparse_sparse_product_impl<Eigen::SparseMatrix<double>,
                                                                 Eigen::SparseMatrix<double>,
                                                                 Eigen::SparseMatrix<double>>(
            mLevelExtractionOps.back(), D, temp, true );

        mLevelExtractionOps.emplace_back( temp.rows() + mActiveFunctions.at( i ).size(), temp.cols() );
        util::verticalConcatInto( temp, active_mat( i ), mLevelExtractionOps.back() );
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

size_t numNonZerosInColumns(const Eigen::SparseMatrix<double>& mat, const std::vector<FunctionId>& cols )
{
    size_t out = 0;
    for( const FunctionId& col : cols )
    {// WARNING: Requires compressed matrix
        out += mat.outerIndexPtr()[col + 1] - mat.outerIndexPtr()[col];
    }
    return out;
}

Eigen::SparseMatrix<double> nonzeroRowsOfColumns( const Eigen::SparseMatrix<double>& mat, const std::vector<FunctionId>& cols )
{
    std::set<Eigen::Index> nonzero_rows;
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve( numNonZerosInColumns( mat, cols ) );

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
    if( c.dim() != mBasisComplex->parametricAtlas().cmap().dim() )
        throw std::invalid_argument( "HierarchicalTPSplineSpace::extractionOperator only works for elements." );

    const auto [ level, level_d ] = mBasisComplex->parametricAtlas().cmap().unrefinedAncestorDartOfCell( c );
    const topology::Cell level_c( level_d, c.dim() );
    const std::vector<FunctionId> level_conn = mRefinementLevels.at( level )->connectivity( level_c );
    const Eigen::MatrixXd level_op = mRefinementLevels.at( level )->extractionOperator( level_c );
    return nonzeroRowsOfColumns( mLevelExtractionOps.at( level ), level_conn ) * level_op;
}

std::vector<FunctionId> HierarchicalTPSplineSpace::connectivity( const topology::Cell& c ) const
{
    const auto [ level, level_d ] = [&]() {
        if( c.dim() == mBasisComplex->parametricAtlas().cmap().dim() )
        {
            return mBasisComplex->parametricAtlas().cmap().unrefinedAncestorDartOfCell( c );
        }
        else
        {
            return mBasisComplex->parametricAtlas().cmap().dartRanges().toLocalDart( c.dart() );
        }
    }();

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
    return mLevelExtractionOps.back().rows();
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

    Eigen::SparseMatrix<double> prolongationOperator( const HierarchicalTPSplineSpace& hier_ss )
    {
        const std::vector<std::vector<FunctionId>>& active_funcs = hier_ss.activeFunctions();
        const auto& refinement_levels = hier_ss.refinementLevels();
        const auto active_mat = [&]( const size_t level_ii ) {
            Eigen::SparseMatrix<double> A( active_funcs.at( level_ii ).size(), refinement_levels.at( level_ii )->numFunctions() );
            A.reserve( Eigen::VectorXi::Constant( A.cols(), 1 ) );
            size_t row_ii = 0;
            for( const FunctionId& fid : active_funcs.at( level_ii ) )
                A.coeffRef( row_ii++, fid ) = 1.0;

            return A;
        };

        Eigen::SparseMatrix<double> D_accum = util::sparseIdentity( refinement_levels.at( 0 )->numFunctions() );
        Eigen::SparseMatrix<double> A_stack = active_mat( 0 );
        for( size_t i = 1; i < refinement_levels.size(); i++ )
        {
            const auto D = refinementOp( *refinement_levels.at( i - 1 ), *refinement_levels.at( i ), 1e-10 );
            D_accum = ( D_accum * D ).eval();
            A_stack = util::verticalConcat( A_stack, active_mat( i ) * D_accum.transpose() );
        }

        return A_stack;
    }

    Eigen::MatrixXd grevillePoints( const HierarchicalTPSplineSpace& ss )
    {
        const auto& refinement_levels = ss.refinementLevels();

        Eigen::MatrixXd out( ss.numFunctions(), ss.basisComplex().parametricAtlas().cmap().dim() );

        std::vector<Eigen::MatrixXd> level_greville_points( refinement_levels.size() );
        for( size_t i = 0; i < refinement_levels.size(); i++ )
        {
            level_greville_points.at( i ) = grevillePoints( *refinement_levels.at( i ) );
        }
        size_t ii = 0;
        size_t jj = 0;
        for( const auto& active_funcs : ss.activeFunctions() )
        {
            for( const auto& f : active_funcs )
            {
                out.row( ii++ ) = level_greville_points.at( jj ).row( f );
            }
            jj++;
        }

        return out;
    }
}