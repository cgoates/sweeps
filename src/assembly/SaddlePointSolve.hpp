#pragma once
#include <Eigen/Sparse>

namespace assembly
{
    int solveMixedEigenproblem( int argc,
                                char** argv,
                                const size_t nev,
                                const Eigen::SparseMatrix<double>& A_in,
                                const Eigen::SparseMatrix<double>& M_in,
                                const Eigen::SparseMatrix<double>& B_in );
}