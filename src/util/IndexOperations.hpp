#pragma once
#include <SmallVector.hpp>
#include <Eigen/Dense>
#include <functional>
#include <variant>

namespace util
{
    using IndexVec = SmallVector<size_t, 3>;
    size_t flatten( const IndexVec& tp_indices, const IndexVec& lengths );

    IndexVec unflatten( const size_t index, const IndexVec& lengths );

    IndexVec dropIndex( const IndexVec& vec, const size_t index_to_drop );

    void iterateTensorProduct( const IndexVec& lengths, const std::function<void( const IndexVec& )>& callback );

    /// @brief Iteration of tensor products in various orders and directions.
    /// @param lengths The number of items in each dimension of the tensor product
    /// @param order The order of fastest increase. For example in 3d, {0,1,2} gives standard tensor product ordering.
    /// {1, 2, 0} would instead increase index 1 the fastest, then index 2, then index 0.
    /// @param direction The direction of iteration in each direction, or a constant value.  True is forward, false is backward.
    /// @param callback A function to call on each tensor product index vector.
    void iterateTensorProduct( const IndexVec& lengths,
                               const util::IndexVec& order,
                               const SmallVector<std::variant<bool, size_t>, 3>& direction,
                               const std::function<void( const IndexVec& )>& callback );

    /// @brief Iteration of tensor product indices in various orders.
    /// See documentation for the version with additional arguments. The direction used here is 'true' in every slot.
    void iterateTensorProduct( const IndexVec& lengths,
                               const util::IndexVec& order,
                               const std::function<void( const IndexVec& )>& callback );

    /// @brief Iterates two compatible tensor products together.
    /// See documentation for iterateTensorProduct above.
    void iterateTensorProductSynchronized( const IndexVec& lengths1,
                                           const IndexVec& lengths2,
                                           const util::IndexVec& order1,
                                           const util::IndexVec& order2,
                                           const SmallVector<std::variant<bool, size_t>, 3>& direction1,
                                           const SmallVector<std::variant<bool, size_t>, 3>& direction2,
                                           const std::function<void( const IndexVec&, const IndexVec& )>& callback );

    void iterateVTKTPOrdering( const IndexVec& lengths, const std::function<void( const size_t& )>& callback );

    Eigen::MatrixXd tensorProduct( const SmallVector<Eigen::VectorXd, 3>& vecs );
}