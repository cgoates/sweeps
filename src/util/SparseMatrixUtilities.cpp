#include <SparseMatrixUtilities.hpp>
#include <iostream>

namespace util
{
    using namespace Eigen;

    // See https://stackoverflow.com/a/41777589
    Eigen::SparseMatrix<double> verticalConcat( const Eigen::SparseMatrix<double>& mat1,
                                                const Eigen::SparseMatrix<double>& mat2 )
    {
        SparseMatrix<double> M( mat1.rows() + mat2.rows(), mat1.cols() );
        verticalConcatInto( mat1, mat2, M );
        return M;
    }

    void verticalConcatInto( const Eigen::SparseMatrix<double>& mat1,
                             const Eigen::SparseMatrix<double>& mat2,
                             Eigen::SparseMatrix<double>& out )
    {
        if( mat1.cols() != mat2.cols() )
        {
            std::cerr << "Matrix 1 has size " << mat1.rows() << ", " << mat1.cols() << ".\n";
            std::cerr << "Matrix 2 has size " << mat2.rows() << ", " << mat2.cols() << ".\n";
            throw std::invalid_argument( "Input matrices must have the same number of columns." );
        }
        if( out.cols() != mat1.cols() or out.rows() != mat1.rows() + mat2.rows() or out.nonZeros() != 0 )
        {
            out = SparseMatrix<double>( mat1.rows() + mat2.rows(), mat1.cols() );
        }


        out.reserve( mat1.nonZeros() + mat2.nonZeros() );
        for( Index c = 0; c < mat1.cols(); ++c )
        {
            out.startVec( c ); // Important: Must be called once for each column before inserting!
            for( SparseMatrix<double>::InnerIterator itL( mat1, c ); itL; ++itL )
                out.insertBack( itL.row(), c ) = itL.value();
            for( SparseMatrix<double>::InnerIterator itC( mat2, c ); itC; ++itC )
                out.insertBack( itC.row() + mat1.rows(), c ) = itC.value();
        }
        out.finalize();
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