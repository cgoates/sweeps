#pragma once
#include <TPSplineSpace.hpp>

namespace fitting
{
    /// @brief Use Coons patch interpolation to create a set of coefficients for a tensor product spline from
    ///        boundary coefficients.
    /// @param ss  The spline space for which to calculate the coefficients.
    /// @param boundary_cpts  The coefficient sets for each of the 2*param_dim boundaries of the manifold.
    /// Typically boundary coefficient sets should agree on shared boundaries, but this is not required; if they don't
    /// agree, lower dimensional boundaries are just defined as the average of the provided values.
    Eigen::MatrixXd coonsPatch( const basis::TPSplineSpace& ss, const SmallVector<Eigen::MatrixXd, 6>& boundary_cpts );
}