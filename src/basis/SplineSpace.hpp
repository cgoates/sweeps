#pragma once
#include <FunctionId.hpp>

namespace topology
{
    class Cell;
}

namespace basis
{
    class BasisComplex;

    class SplineSpace
    {
        public:
        virtual const BasisComplex& basisComplex() const = 0;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const = 0;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const = 0;

        virtual size_t numFunctions() const = 0;
    };
}