#include <AABB.hpp>
#include <SimplexUtilities.hpp>

namespace mapping
{
    std::ostream& operator<<( std::ostream& os, const AABB& aabb )
    {
        os << "AABB(" << aabb.min().transpose() << "; " << aabb.max().transpose() << ")";
        return os;
    }

    template<int DIM>
    AABB aabbFromTriangle( const Triangle<DIM>& tri )
    {
        const Vector3dMax mins = tri.v1.array().min( tri.v2.array() ).min( tri.v3.array() );
        const Vector3dMax maxs = tri.v1.array().max( tri.v2.array() ).max( tri.v3.array() );
        return AABB( mins, maxs );
    }
    template AABB aabbFromTriangle<2>( const Triangle<2>& tri );
    template AABB aabbFromTriangle<3>( const Triangle<3>& tri );
}