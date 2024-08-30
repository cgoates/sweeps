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

        const SmallVector<TPSplineSpace, 3>& scalarTPBases() const { return mScalarTPBases; }

        private:
        const std::shared_ptr<const DivConfBasisComplex> mBasisComplex;

        // These are stored just for ownership and lifetime purposes.
        SmallVector<std::shared_ptr<const BasisComplex1d>, 3> mReducedDegree1dBasisComplex;
        SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> mReducedDegree1dBases;
        SmallVector<std::shared_ptr<const TPBasisComplex>, 3> mScalarTPBasisComplexes;
        SmallVector<std::shared_ptr<const TPSplineSpace>, 3> m2dSourceTPBases; // only used in 3d

        // These are the scalar bases we use directly to create the vector basis.
        SmallVector<TPSplineSpace, 3> mScalarTPBases;
    };
}