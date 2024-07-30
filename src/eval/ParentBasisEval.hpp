#pragma once
#include <Eigen/Dense>

namespace param
{
    class ParentPoint;
}

namespace basis
{
    class ParentBasis;
}

namespace eval
{
    class ParentBasisEval
    {
        public:
        ParentBasisEval( const basis::ParentBasis& pb, const param::ParentPoint& point, const size_t n_derivatives );
        ParentBasisEval() {}

        Eigen::MatrixXd mEvals;
    };
}