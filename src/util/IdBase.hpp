#pragma once
#include <cstddef>
#include <ostream>

namespace Eigen
{
    using Index = std::ptrdiff_t;
}

namespace util
{
    /// @brief An id class that provides type-safe ids.  See https://stackoverflow.com/a/63821945
    /// @tparam T  A tag that differentiates ids from each other.
    /// As a usage example,
    /// struct FunctionIdTag;
    /// using FunctionId = IdBase<FunctionIdTag>;
    template<typename T>
    class IdBase
    {
        public:
        IdBase( const Eigen::Index id ) : mId( id ) {}
        operator Eigen::Index() const { return mId; }

        Eigen::Index id() const { return mId; }

        bool operator<( const IdBase& o ) const
        {
            return mId < o.id();
        }

        bool operator>( const IdBase& o ) const
        {
            return mId > o.id();
        }

        friend std::ostream& operator<<( std::ostream& os, const IdBase& id )
        {
            os << id.mId;
            return os;
        }

        private:
        Eigen::Index mId;
    };
}