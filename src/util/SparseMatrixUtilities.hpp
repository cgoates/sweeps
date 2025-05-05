#pragma once
#include <Eigen/Sparse>

namespace util
{
    Eigen::SparseMatrix<double> verticalConcat( const Eigen::SparseMatrix<double>& mat1,
                                                const Eigen::SparseMatrix<double>& mat2 );

    // Concatenates two matrices vertically into a provided output matrix.
    // Can be more performant than the version that returns a matrix in some scenarios.
    void verticalConcatInto( const Eigen::SparseMatrix<double>& mat1,
                             const Eigen::SparseMatrix<double>& mat2,
                             Eigen::SparseMatrix<double>& out );

    Eigen::SparseMatrix<double> sparseIdentity( const Eigen::Index n );
}