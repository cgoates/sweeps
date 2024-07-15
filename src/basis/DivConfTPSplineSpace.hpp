#pragma once
#include <SplineSpace.hpp>
#include <DivConfBasisComplex.hpp>
#include <TPSplineSpace.hpp>

namespace basis
{
    class DivConfTPSplineSpace : public SplineSpace
    {
        public:
        DivConfTPSplineSpace( const DivConfBasisComplex& bc, const TPSplineSpace& primal_basis );

        virtual const DivConfBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        private:
        const DivConfBasisComplex& mBasisComplex;

        // These are stored just for ownership and lifetime purposes.
        std::vector<BasisComplex1d> mReducedDegree1dBasisComplex;
        std::vector<BSplineSpace1d> mReducedDegree1dBases;
        std::vector<TPBasisComplex> mScalarTPBasisComplexes;
        std::vector<TPSplineSpace> m2dSourceTPBases; // only used in 3d

        // These are the scalar bases we use directly to create the vector basis.
        std::vector<TPSplineSpace> mScalarTPBases;
    };
}