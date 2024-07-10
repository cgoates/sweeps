#pragma once
#include <Eigen/Dense>

namespace basis
{
    class FunctionId
    {
        public:
        explicit FunctionId( const Eigen::Index id ) : mId( id ) {}
        operator Eigen::Index() const { return mId; }

        Eigen::Index id() const { return mId; }

        bool operator<( const FunctionId& o ) const
        {
            return mId < o.id();
        }

        bool operator>( const FunctionId& o ) const
        {
            return mId > o.id();
        }

        bool operator==( const FunctionId& o ) const
        {
            return mId == o.id();
        }

        friend std::ostream& operator<<( std::ostream& os, const FunctionId& fid )
        {
            os << fid.mId;
            return os;
        }

        private:
        Eigen::Index mId;
    };
}