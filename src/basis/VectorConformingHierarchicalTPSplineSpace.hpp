#pragma once
#include <SplineSpace.hpp>
#include <VectorConformingBasisComplex.hpp>
#include <HierarchicalTPSplineSpace.hpp>

namespace basis
{
    class VectorConformingHierarchicalTPSplineSpace : public SplineSpace
    {
        public:
        VectorConformingHierarchicalTPSplineSpace( const std::shared_ptr<const VectorConformingBasisComplex>& bc, const HierarchicalTPSplineSpace& primal_basis );

        virtual const VectorConformingBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        virtual size_t numVectorComponents() const override { return mScalarTPBases.size(); }

        const SmallVector<std::shared_ptr<const HierarchicalTPSplineSpace>, 3>& scalarBases() const { return mScalarTPBases; }

        private:
        const std::shared_ptr<const VectorConformingBasisComplex> mBasisComplex;
        SmallVector<std::shared_ptr<const HierarchicalTPSplineSpace>, 3> mScalarTPBases;
    };
}