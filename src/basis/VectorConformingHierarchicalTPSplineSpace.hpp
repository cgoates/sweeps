#pragma once
#include <SplineSpace.hpp>
#include <VectorConformingBasisComplex.hpp>
#include <VectorConformingTPSplineSpace.hpp>
#include <HierarchicalTPSplineSpace.hpp>

namespace basis
{
    class VectorConformingHierarchicalTPSplineSpace : public SplineSpace
    {
        public:
        VectorConformingHierarchicalTPSplineSpace( const std::shared_ptr<const VectorConformingBasisComplex>& bc, const HierarchicalTPSplineSpace& primal_basis );
        VectorConformingHierarchicalTPSplineSpace(
            const std::shared_ptr<const VectorConformingBasisComplex>& bc,
            const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>& refinement_levels,
            const std::vector<std::vector<FunctionId>>& active_funcs );
        virtual ~VectorConformingHierarchicalTPSplineSpace() = default;

        virtual const VectorConformingBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        virtual size_t numVectorComponents() const override { return mScalarTPBases.size(); }

        /// Local hierarchical vector functions are ordered component-major, level-minor:
        /// all active functions for component 0 across every refinement level, then
        /// component 1, and so on. This matches connectivity() and must be preserved by
        /// patch-local to global maps in hierarchical multipatch vector spaces.
        const SmallVector<std::shared_ptr<const HierarchicalTPSplineSpace>, 3>& scalarBases() const { return mScalarTPBases; }

        private:
        const std::shared_ptr<const VectorConformingBasisComplex> mBasisComplex;
        SmallVector<std::shared_ptr<const HierarchicalTPSplineSpace>, 3> mScalarTPBases;
    };
}
