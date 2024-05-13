#pragma once
#include <Eigen/Dense>

namespace basis
{
    class FunctionId
    {
        public:
        Eigen::Index operator()() const { return mId; }

        bool operator<( const FunctionId& o ) const
        {
            return mId < o();
        }

        bool operator>( const FunctionId& o ) const
        {
            return mId > o();
        }

        bool operator==( const FunctionId& o ) const
        {
            return mId == o();
        }

        private:
        Eigen::Index mId;
    };
}