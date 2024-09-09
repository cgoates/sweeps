#pragma once
#include <ParentDomain.hpp>
#include <CustomEigen.hpp>

namespace param
{
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
    ParentPoint compressCoordinates( const ParentDomain& domain, const Vector6dMax& coords, const BaryCoordIsZeroVec& zeros );
    ParentPoint compressCoordinates( const ParentDomain& domain, const Vector6dMax& expanded_coords, const double is_zero_tol );
    ParentPoint average( const ParentPoint& pt1, const ParentPoint& pt2 );
    ParentPoint tensorProduct( const ParentPoint& pt1, const ParentPoint& pt2 );
    ParentPoint liftFromBoundary( const ParentPoint& bdry_point, const ParentDomain& interior_domain, const BaryCoordIsZeroVec& bdry );
    std::ostream& operator<<( std::ostream& o, const ParentPoint& ppt );
    BaryCoordIsZeroVec join( const BaryCoordIsZeroVec& v1, const BaryCoordIsZeroVec& v2 );
}