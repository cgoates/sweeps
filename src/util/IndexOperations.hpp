#pragma once
#include <SmallVector.hpp>
#include <Eigen/Dense>
#include <functional>

namespace util
{
    using IndexVec = SmallVector<size_t, 3>;
    size_t flatten( const IndexVec& tp_indices, const IndexVec& lengths );

    IndexVec unflatten( const size_t index, const IndexVec& lengths );

    void iterateTensorProduct( const IndexVec& lengths, const std::function<void( const IndexVec& )>& callback );

    void iterateVTKTPOrdering( const IndexVec& lengths, const std::function<void( const size_t& )>& callback );

    Eigen::MatrixXd tensorProduct( const SmallVector<Eigen::VectorXd, 3>& vecs );
}