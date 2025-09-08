#pragma once
#include <Eigen/Sparse>
#include <slepceps.h>
#include <vector>
#include <tuple>

namespace slepc_utils
{
    PetscErrorCode eigenSparseToPetsc( const Eigen::SparseMatrix<double>& eigenMat,
                                       Mat& petscMat,
                                       MPI_Comm comm = PETSC_COMM_WORLD );

    // Solve generalized eigenvalue problem and return results
    std::pair<std::vector<PetscScalar>, std::vector<Eigen::VectorXd>>
        solveGeneralizedEigenvalueSparse( const Eigen::SparseMatrix<double>& K1_bar,
                                          const Eigen::SparseMatrix<double>& Q,
                                          const Eigen::SparseMatrix<double>& K2,
                                          int nev,
                                          MPI_Comm comm = PETSC_COMM_WORLD );

    std::pair<std::vector<PetscScalar>, std::vector<Eigen::VectorXd>>
        solveGeneralizedEigenvalueSparse( const Eigen::SparseMatrix<double>& A,
                                          const Eigen::SparseMatrix<double>& B,
                                          int nev,
                                          int n_edge,
                                          MPI_Comm comm = PETSC_COMM_WORLD );
} // namespace slepc_utils