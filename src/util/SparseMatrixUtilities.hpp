#pragma once
#include <Eigen/Sparse>

namespace util
{
    Eigen::SparseMatrix<double> verticalConcat( const Eigen::SparseMatrix<double>& mat1,
                                                const Eigen::SparseMatrix<double>& mat2 );

    Eigen::SparseMatrix<double> sparseIdentity( const Eigen::Index n );
}