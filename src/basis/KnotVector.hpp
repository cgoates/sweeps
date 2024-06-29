#pragma once
#include <vector>
#include <Eigen/Dense>

namespace basis
{
    class KnotVector
    {
        public:
        KnotVector( const std::vector<double>& knots, const double parametric_tol );
        size_t size() const;
        double knot( const size_t knot_ii ) const;
        std::vector<double> uniqueKnots() const;
        private:
        std::vector<std::pair<double, size_t>> mKnots;
    };

    Eigen::VectorXd parametricLengths( const KnotVector& kv );
}