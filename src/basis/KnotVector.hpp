#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SmallVector.hpp>

namespace basis
{
    class KnotVector
    {
        public:
        KnotVector( const std::vector<double>& knots, const double parametric_tol );
        KnotVector( const std::vector<std::pair<double, size_t>>& knot_multiplicities );
        size_t size() const;
        double knot( const size_t knot_ii ) const;
        std::vector<double> uniqueKnots() const;

        void iterateUniqueKnots( const std::function<void( const double, const size_t )>& iter ) const;
        const std::vector<std::pair<double, size_t>>& uniqueKnotMultiplicities() const { return mKnots; }

        void duplicate( const size_t knot_ii );
        void insert( const double knot, const double param_tol );

        private:
        std::vector<std::pair<double, size_t>> mKnots;
    };

    std::ostream& operator<<( std::ostream& o, const KnotVector& kv );

    /// @brief Gives the parametric lengths of the cells represented by the knot vector.
    /// @param kv The knot vector
    /// @return The differences between each unique knot value
    Eigen::VectorXd parametricLengths( const KnotVector& kv );

    /// @brief Gives the number of parametric elements in the knot vector with non-zero width.
    /// @param kv The knot vector
    /// @return Number of nonzero-width knot intervals
    size_t numElements( const KnotVector& kv );

    /// @brief Gives the extraction operator between the given knot vector and one that has all knots duplicated to
    /// multiplicity >= degree
    /// @param kv  The knot vector
    /// @param degree  The degree of the spline represented by the knot vector
    /// @return  The extraction operator
    Eigen::SparseMatrix<double> globalExtractionOp( const KnotVector& kv, const size_t degree );

    /// @brief Gives an expression of the coarse basis in terms of the fine basis.
    /// @param kv_coarse The knot vector representing the coarse basis.
    /// @param kv_fine The knot vector representing the fine basis.  Must be nested within kv_coarse.
    /// @param degree The degree of both bases.
    /// @return An extraction operator of the coarse basis in terms of the fine basis.
    Eigen::SparseMatrix<double> refinementOp( const KnotVector& kv_coarse,
                                              const KnotVector& kv_fine,
                                              const size_t degree,
                                              const double param_tol );

    /// @brief Gives an expression of the coarse basis in terms of the fine basis for a tensor product B-spline.
    /// @param kvs_coarse The knot vectors representing the coarse basis.
    /// @param kvs_fine The knot vectors representing the fine basis.  Must be nested within kvs_coarse.
    /// @param degrees The degrees of both bases.
    /// @return An extraction operator of the coarse basis in terms of the fine basis.
    Eigen::SparseMatrix<double> refinementOp( const SmallVector<KnotVector, 3>& kvs_coarse,
                                              const SmallVector<KnotVector, 3>& kvs_fine,
                                              const SmallVector<size_t, 3> degrees,
                                              const double param_tol );

    /// @brief Reduces the order of the knot vector, essentially removing one of the beginning and end knots.
    /// @param kv The knot vector to reduce the order of
    /// @return The reduced order knot vector
    KnotVector reducedOrder( const KnotVector& kv );

    /// @brief Calculate the greville points from a knot vector
    /// @param kv  The knot vector
    /// @param degree  The basis degree
    /// @return  The greville points
    Eigen::VectorXd grevillePoints( const KnotVector& kv, const size_t degree );

    KnotVector integerKnotsWithNElems( const size_t n_elems, const size_t degree );

    KnotVector dyadicRefine( const KnotVector& kv );
    KnotVector nAdicRefine( const KnotVector& kv, const size_t n );
}