#pragma once
#include <SplineSpace.hpp>
#include <DivConfBasisComplex.hpp>
#include <TPSplineSpace.hpp>

namespace basis
{
    class DivConfTPSplineSpace : public SplineSpace
    {
        public:
        DivConfTPSplineSpace( const std::shared_ptr<const DivConfBasisComplex>& bc, const TPSplineSpace& primal_basis );

        virtual const DivConfBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        virtual size_t numVectorComponents() const override { return 2; }

        const SmallVector<std::shared_ptr<const BSplineSpace1d>, 3>& reducedDegree1dBases() const
        {
            return mReducedDegree1dBases;
        }

        const SmallVector<std::shared_ptr<const TPSplineSpace>, 3>& scalarTPBases() const { return mScalarTPBases; }

        private:
        const std::shared_ptr<const DivConfBasisComplex> mBasisComplex;
        SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> mReducedDegree1dBases;
        SmallVector<std::shared_ptr<const TPSplineSpace>, 3> mScalarTPBases;
    };
}