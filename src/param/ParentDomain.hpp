#pragma once
#include <SmallVector.hpp>
#include <Eigen/Dense>
#include <numeric>

using Vector3dMax = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3>;
using Vector6dMax = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 6>;

namespace param
{
    /// @brief Represents a barycentric coordinate system of a given dimension.
    class CoordinateSystem
    {
        public:
        CoordinateSystem() {}
        CoordinateSystem( const size_t dim ) : mDim( dim ) {}
        size_t dim() const { return mDim; }

        private:
        size_t mDim;
    };

    size_t numTotalCoordinates( const CoordinateSystem& cs );

    /// @brief Represents a parent domain coordinate system made up of tensor products
    /// of one or more barycentric coordinate systems.
    class ParentDomain
    {
        public:
        ParentDomain( std::span<const CoordinateSystem> cs_vec ) : mCoordinateGroups( cs_vec ) {}
        const SmallVector<CoordinateSystem, 3>& coordinateGroups() const { return mCoordinateGroups; }
        private:
        SmallVector<CoordinateSystem, 3> mCoordinateGroups;
    };

    ParentDomain simplexDomain( const size_t dim );
    ParentDomain cubeDomain( const size_t dim );

    size_t numGroups( const ParentDomain& pd );
    size_t dim( const ParentDomain& pd );
    size_t numTotalCoordinates( const ParentDomain& pd );
    void iterateGroups( const ParentDomain& pd, const std::function<void( const size_t, const size_t, const CoordinateSystem& )>& callback );
    Vector6dMax expandedCoordinates( const ParentDomain& domain, const Vector3dMax& pt );

    using BaryCoordIsZeroVec = SmallVector<bool, 6>;

    /// @brief Represents a point in a parent domain
    class ParentPoint
    {
        public:
        ParentPoint( const ParentDomain& domain, const Vector3dMax& point, const BaryCoordIsZeroVec& zero_vec );
        ParentDomain mDomain;
        Vector3dMax mPoint;
        BaryCoordIsZeroVec mBaryCoordIsZero;
    };

    ParentPoint pointOnBoundary( const ParentDomain& domain, const BaryCoordIsZeroVec& is_zero );
    Vector6dMax expandedCoordinates( const ParentPoint& pt );
    ParentPoint compressCoordinates( const ParentDomain& domain, const Vector6dMax& expanded_coords, const double is_zero_tol );
    BaryCoordIsZeroVec join( const BaryCoordIsZeroVec& v1, const BaryCoordIsZeroVec& v2 );
}