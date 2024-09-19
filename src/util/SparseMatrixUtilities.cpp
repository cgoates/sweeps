#include <SparseMatrixUtilities.hpp>
#include <iostream>

namespace util
{
    Eigen::SparseMatrix<double> verticalConcat( const Eigen::SparseMatrix<double>& mat1,
                                                const Eigen::SparseMatrix<double>& mat2 )
    {
        if( mat1.cols() != mat2.cols() )
        {
            std::cerr << "Matrix 1 has size " << mat1.rows() << ", " << mat1.cols() << ".\n";
            std::cerr << "Matrix 2 has size " << mat2.rows() << ", " << mat2.cols() << ".\n";
            throw std::invalid_argument( "Input matrices must have the same number of columns." );
        }

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve( mat1.nonZeros() + mat2.nonZeros() );

        for( Eigen::Index k = 0; k < mat1.outerSize(); ++k )
        {
            for( Eigen::SparseMatrix<double>::InnerIterator it( mat1, k ); it; ++it )
            {
                triplets.emplace_back( it.row(), it.col(), it.value() );
            }
        }

        for( Eigen::Index k = 0; k < mat2.outerSize(); ++k )
        {
            for( Eigen::SparseMatrix<double>::InnerIterator it( mat2, k ); it; ++it )
            {
                triplets.emplace_back( it.row() + mat1.rows(), it.col(), it.value() );
            }
        }

        Eigen::SparseMatrix<double> result( mat1.rows() + mat2.rows(), mat1.cols() );
        result.setFromTriplets( triplets.begin(), triplets.end() );

        return result;
    }

    Eigen::SparseMatrix<double> sparseIdentity( const Eigen::Index n )
    {
        Eigen::SparseMatrix<double> out( n, n );
        out.reserve( Eigen::VectorXi::Ones( n ) );

        for( Eigen::Index i = 0; i < n; i++ )
            out.coeffRef( i, i ) = 1;

        return out;
    }
}