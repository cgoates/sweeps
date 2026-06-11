#pragma once
#include <FunctionId.hpp>
#include <KnotVector.hpp>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

namespace basis
{
    struct LocalExtraction
    {
        std::vector<FunctionId> connectivity;
        Eigen::MatrixXd extraction;
    };

    /// @brief Gives element-local extraction operators for a univariate B-spline space.
    /// @param kv The knot vector defining the space.
    /// @param degree The B-spline degree.
    /// @return One local extraction per nonzero-width knot interval. Rows are local B-spline functions and columns are
    /// local Bernstein functions on the interval.
    std::vector<LocalExtraction> localExtractions( const KnotVector& kv, const size_t degree );

    /// @brief Gives the starting column in globalExtractionOp(kv, degree) for each element.
    ///
    /// globalExtractionOp inserts every interior knot until its multiplicity is at least degree, producing a C0-style
    /// Bezier column layout for all spaces with interior multiplicity <= degree. Fully discontinuous C^{-1} spaces
    /// already have interior multiplicity degree + 1, so adjacent elements do not share an endpoint column.
    std::vector<size_t> elementBezierColumnOffsets( const KnotVector& kv, const size_t degree );

    /// @brief Gives the Bernstein degree elevation matrix from source_degree to target_degree.
    /// @return Matrix M such that B_source = M * B_target.
    Eigen::MatrixXd bernsteinDegreeElevationMatrix( const size_t source_degree, const size_t target_degree );

    /// @brief Gives the knot vector for degree elevation by increasing every knot multiplicity.
    /// @param kv The source knot vector.
    /// @param elevation_amount The degree increase.
    /// @return The elevated knot vector.
    KnotVector degreeElevatedKnotVector( const KnotVector& kv, const size_t elevation_amount );

    /// @brief Gives the exact transfer operator for univariate B-spline degree elevation.
    /// @return Operator T such that elevated_control_points = source_control_points * T.
    Eigen::SparseMatrix<double> degreeElevationOp( const KnotVector& kv,
                                                   const size_t source_degree,
                                                   const size_t target_degree,
                                                   const double param_tol );
}
