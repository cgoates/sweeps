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

        const std::vector<std::shared_ptr<const BSplineSpace1d>>& reducedDegree1dBases() const
        {
            return mReducedDegree1dBases;
        }

        private:
        const std::shared_ptr<const DivConfBasisComplex> mBasisComplex;

        // These are stored just for ownership and lifetime purposes.
        std::vector<std::shared_ptr<const BasisComplex1d>> mReducedDegree1dBasisComplex;
        std::vector<std::shared_ptr<const BSplineSpace1d>> mReducedDegree1dBases;
        std::vector<std::shared_ptr<const TPBasisComplex>> mScalarTPBasisComplexes;
        std::vector<std::shared_ptr<const TPSplineSpace>> m2dSourceTPBases; // only used in 3d

        // These are the scalar bases we use directly to create the vector basis.
        std::vector<TPSplineSpace> mScalarTPBases;
    };
}