#pragma once
#include <MultiPatchSplineSpace.hpp>
#include <HierarchicalMultiPatchBasisComplex.hpp>
#include <HierarchicalTPSplineSpace.hpp>

namespace basis
{
    class HierarchicalMultiPatchSplineSpace : public SplineSpace
    {
        public:
        HierarchicalMultiPatchSplineSpace( const std::shared_ptr<const HierarchicalMultiPatchBasisComplex>& bc,
                                   const std::vector<std::shared_ptr<const MultiPatchSplineSpace>>& refinement_levels );
        virtual ~HierarchicalMultiPatchSplineSpace() = default;

        virtual const HierarchicalMultiPatchBasisComplex& basisComplex() const override;
        const std::shared_ptr<const HierarchicalMultiPatchBasisComplex>& basisComplexPtr() const;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        const std::vector<std::shared_ptr<const MultiPatchSplineSpace>>& refinementLevels() const { return mRefinementLevels; }
        const std::vector<std::shared_ptr<const HierarchicalTPSplineSpace>>& constituents() const { return mConstituents; }
        const std::vector<std::vector<FunctionId>>& functionIdMap() const { return mFuncIds; }

        private:

        const std::shared_ptr<const HierarchicalMultiPatchBasisComplex> mBasisComplex;
        const std::vector<std::shared_ptr<const MultiPatchSplineSpace>> mRefinementLevels;
        std::vector<std::shared_ptr<const HierarchicalTPSplineSpace>> mConstituents;
        std::vector<std::vector<FunctionId>> mFuncIds; // ith constituent's jth function has id mFuncIds[i][j] in the multipatch space
        size_t mNumActiveFuncs;
    };

    HierarchicalMultiPatchSplineSpace
        buildHierarchicalSplineSpace( const std::vector<std::shared_ptr<const MultiPatchSplineSpace>>& refinement_levels,
                                      const std::vector<std::vector<topology::Cell>>& leaf_elements );
}