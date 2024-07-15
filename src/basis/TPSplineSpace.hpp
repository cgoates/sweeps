#pragma once
#include <BSplineSpace1d.hpp>
#include <TPBasisComplex.hpp>

namespace basis
{
    class TPSplineSpace : public SplineSpace
    {
        public:
        TPSplineSpace( const TPBasisComplex& bc, const SplineSpace& source, const BSplineSpace1d& line );

        virtual const TPBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        const BSplineSpace1d& line() const { return mLine; }
        const SplineSpace& source() const { return mSource; }

        private:

        FunctionId flatten( const FunctionId& source_fid, const FunctionId& line_fid ) const;

        const TPBasisComplex& mBasisComplex;
        const SplineSpace& mSource;
        const BSplineSpace1d& mLine;
    };
}