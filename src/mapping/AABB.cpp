#include <AABB.hpp>

namespace mapping
{
    std::ostream& operator<<( std::ostream& os, const AABB& aabb )
    {
        os << "AABB(" << aabb.min().transpose() << "; " << aabb.max().transpose() << ")";
        return os;
    }
}