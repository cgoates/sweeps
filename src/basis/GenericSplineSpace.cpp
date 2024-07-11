#include <GenericSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParametricAtlas.hpp>
#include <BasisComplex1d.hpp>
#include <Eigen/Sparse>
#include <KnotVector.hpp>

namespace basis
{
    GenericSplineSpace::GenericSplineSpace( const BasisComplex& bc,
                                            const std::map<size_t, Eigen::MatrixXd>& ex_ops,
                                            const std::map<size_t, std::vector<FunctionId>>& connectivity )
        : mBasisComplex( bc ), mExtractionOps( ex_ops ), mConnectivity( connectivity )
    {
        mNumFunctions = 0;
        for( const auto& pr : connectivity )
        {
            for( const FunctionId& fid : pr.second ) mNumFunctions = std::max( mNumFunctions, size_t( fid ) + 1 );
        }
    }

    const BasisComplex& GenericSplineSpace::basisComplex() const
    {
        return mBasisComplex;
    }

    Eigen::MatrixXd GenericSplineSpace::extractionOperator( const topology::Cell& c ) const
    {
        if( c.dim() != mBasisComplex.parametricAtlas().cmap().dim() )
            throw std::runtime_error( "Extraction operators only supported for elements" );
        const auto indexing = indexingOrError( mBasisComplex.parametricAtlas().cmap(), c.dim() );
        return mExtractionOps.at( indexing( c ) );
    }

    std::vector<FunctionId> GenericSplineSpace::connectivity( const topology::Cell& c ) const
    {
        if( c.dim() != mBasisComplex.parametricAtlas().cmap().dim() )
            throw std::runtime_error( "Extraction operators only supported for elements" );
        const auto indexing = indexingOrError( mBasisComplex.parametricAtlas().cmap(), c.dim() );
        return mConnectivity.at( indexing( c ) );
    }

    size_t GenericSplineSpace::numFunctions() const
    {
        return mNumFunctions;
    }


    GenericSplineSpace knotVectorSplineSpace( const BasisComplex1d& bc, const KnotVector& kv )
    {
        const size_t degree = bc.defaultParentBasis().mBasisGroups.at(0).degrees.at(0);
        using SparseMatrixXd = Eigen::SparseMatrix<double>;

        const SparseMatrixXd C = globalExtractionOp( kv, degree );

        const topology::CombinatorialMap1d& cmap = bc.parametricAtlas().cmap();

        const auto elem_ids = indexingOrError( cmap, 1 );
        std::map<size_t, Eigen::MatrixXd> ex_ops;
        std::map<size_t, std::vector<FunctionId>> connectivity;

        const auto get_elem_connectivity = [&degree,&C]( const size_t elem_ii ) {
            std::set<FunctionId> unique_rows;
            for( size_t col_ii = 0; col_ii <= degree; col_ii++ )
            {
                for( SparseMatrixXd::InnerIterator it( C, col_ii + elem_ii * degree ); it; ++it )
                {
                    unique_rows.insert( FunctionId( it.row() ) );
                }
            }
            return std::vector<FunctionId>( unique_rows.begin(), unique_rows.end() );
        };

        const auto get_elem_operator = [&degree, &C]( const size_t elem_ii, const std::vector<FunctionId>& conn_elem ) {
            Eigen::MatrixXd C_elem = Eigen::MatrixXd::Zero( degree + 1, degree + 1 );
            for( size_t col = 0; col <= degree; col++ )
            {
                for( SparseMatrixXd::InnerIterator it( C, col + elem_ii * degree ); it; ++it )
                {
                    const Eigen::Index row = std::distance( conn_elem.begin(), std::find( conn_elem.begin(), conn_elem.end(), it.row() ) );
                    C_elem( row, col ) = it.value();
                }
            }
            return C_elem;
        };

        iterateCellsWhile( cmap, 1, [&]( const topology::Edge& e ) {
            const size_t elem_ii = elem_ids( e );
            const std::vector<FunctionId>& conn_elem =
                connectivity.emplace( elem_ii, get_elem_connectivity( elem_ii ) ).first->second;
            ex_ops.emplace( elem_ii, get_elem_operator( elem_ii, conn_elem ) );
            return true;
        } );

        return GenericSplineSpace( bc, ex_ops, connectivity );
    }
}