#pragma once
#include <FunctionId.hpp>

namespace topology
{
    class Cell;
}

namespace basis
{
    class BasisMesh;

    class SplineSpace
    {
        public:
        virtual const BasisMesh& basisMesh() const = 0;

        virtual const Eigen::MatrixXd extractionOperator( const topology::Cell& ) const = 0;

        virtual const std::vector<FunctionId> connectivity( const topology::Cell& ) const = 0;
    };
}