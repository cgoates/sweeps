#include <SparseMatrixUtilities.hpp>
#include <iostream>

namespace util
{
    // See https://stackoverflow.com/a/41777589
    Eigen::SparseMatrix<double> verticalConcat( const Eigen::SparseMatrix<double>& mat1,
                                                const Eigen::SparseMatrix<double>& mat2 )
    {
        if( mat1.cols() != mat2.cols() )
        {
            std::cerr << "Matrix 1 has size " << mat1.rows() << ", " << mat1.cols() << ".\n";
            std::cerr << "Matrix 2 has size " << mat2.rows() << ", " << mat2.cols() << ".\n";
            throw std::invalid_argument( "Input matrices must have the same number of columns." );
        }
        using namespace Eigen;
        SparseMatrix<double> M( mat1.rows() + mat2.rows(), mat1.cols() );
        M.reserve( mat1.nonZeros() + mat2.nonZeros() );
        for( Index c = 0; c < mat1.cols(); ++c )
        {
            M.startVec( c ); // Important: Must be called once for each column before inserting!
            for( SparseMatrix<double>::InnerIterator itL( mat1, c ); itL; ++itL )
                M.insertBack( itL.row(), c ) = itL.value();
            for( SparseMatrix<double>::InnerIterator itC( mat2, c ); itC; ++itC )
                M.insertBack( itC.row() + mat1.rows(), c ) = itC.value();
        }
        M.finalize();
        return M;
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