#pragma once
#include <SmallVector.hpp>
#include <Eigen/Dense>
#include <numeric>
#include <CustomEigen.hpp>

namespace param
{
    using BaryCoordIsZeroVec = SmallVector<bool, 6>;

    /// @brief Represents a barycentric coordinate system of a given dimension.
    class CoordinateSystem
    {
        public:
        CoordinateSystem() {}
        CoordinateSystem( const size_t dim ) : mDim( dim ) {}
        size_t dim() const { return mDim; }

        bool operator!=( const CoordinateSystem& o ) const
        {
            return dim() != o.dim();
        }

        bool operator==( const CoordinateSystem& o ) const
        {
            return dim() == o.dim();
        }

        private:
        size_t mDim;
    };

    size_t numTotalCoordinates( const CoordinateSystem& cs );
    size_t numParametricLengths( const CoordinateSystem& cs );

    /// @brief Represents a parent domain coordinate system made up of tensor products
    /// of one or more barycentric coordinate systems.
    class ParentDomain
    {
        public:
        ParentDomain( std::span<const CoordinateSystem> cs_vec ) : mCoordinateGroups( cs_vec.begin(), cs_vec.end() ) {}
        const SmallVector<CoordinateSystem, 3>& coordinateGroups() const { return mCoordinateGroups; }

        bool operator!=( const ParentDomain& o ) const
        {
            return coordinateGroups() != o.coordinateGroups();
        }

        bool operator==( const ParentDomain& o ) const
        {
            return coordinateGroups() == o.coordinateGroups();
        }

        private:
        SmallVector<CoordinateSystem, 3> mCoordinateGroups;
    };

    ParentDomain simplexDomain( const size_t dim );
    ParentDomain cubeDomain( const size_t dim );
    ParentDomain tensorProduct( const ParentDomain& pd1, const ParentDomain& pd2 );

    std::ostream& operator<<( std::ostream& o, const ParentDomain& pd );

    size_t numGroups( const ParentDomain& pd );
    size_t dim( const ParentDomain& pd );
    size_t numTotalCoordinates( const ParentDomain& pd );
    size_t numParametricLengths( const ParentDomain& pd );
    bool isCartesian( const ParentDomain& pd );
    /// @brief Iterate the groups of the parent domain.
    /// @param pd Parent domain to iterate.
    /// @param callback Calls back on each group with the first expanded coordinate and explicit coordinate belonging to the group, in that order.
    void iterateGroups( const ParentDomain& pd, const std::function<void( const size_t, const size_t, const CoordinateSystem& )>& callback );
    Vector6dMax expandedCoordinates( const ParentDomain& domain, const Vector3dMax& pt );
    SmallVector<size_t, 6> changingCoordinates( const ParentDomain& pd, const BaryCoordIsZeroVec& bdry );
}