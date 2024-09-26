#pragma once
#include <SplineSpace.hpp>
#include <DivConfBasisComplex.hpp>
#include <HierarchicalTPSplineSpace.hpp>

namespace basis
{
    class DivConfHierarchicalTPSplineSpace : public SplineSpace
    {
        public:
        DivConfHierarchicalTPSplineSpace( const std::shared_ptr<const DivConfBasisComplex>& bc, const HierarchicalTPSplineSpace& primal_basis );

        virtual const DivConfBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        virtual size_t numVectorComponents() const override { return 2; }

        const SmallVector<std::shared_ptr<const HierarchicalTPSplineSpace>, 3>& scalarBases() const { return mScalarTPBases; }

        private:
        const std::shared_ptr<const DivConfBasisComplex> mBasisComplex;
        SmallVector<std::shared_ptr<const HierarchicalTPSplineSpace>, 3> mScalarTPBases;
    };
}