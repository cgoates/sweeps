#pragma once
#include <CustomEigen.hpp>
#include <iosfwd>

namespace mapping
{
    class AABB
    {
        public:
        AABB( const Vector3dMax& min, const Vector3dMax& max ) : mMin( min ), mMax( max ) {}

        const Vector3dMax& min() const { return mMin; }
        const Vector3dMax& max() const { return mMax; }

        bool contains( const Vector3dMax& pt ) const
        {
            return ( pt.array() >= mMin.array() ).all() and ( pt.array() <= mMax.array() ).all();
        }

        private:
        Vector3dMax mMin;
        Vector3dMax mMax;
    };

    std::ostream& operator<<( std::ostream& os, const AABB& aabb );
}