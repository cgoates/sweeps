#pragma once
#include <HierarchicalMultiPatchSplineSpace.hpp>
#include <VectorConformingMultiPatchSplineSpace.hpp>
#include <VectorConformingHierarchicalTPSplineSpace.hpp>

namespace basis
{
    struct VectorHierarchicalLevelFunctionId
    {
        size_t level;
        FunctionId function;
        bool orientation;
    };

    class VectorConformingHierarchicalMultiPatchSplineSpace : public SplineSpace
    {
        public:
        VectorConformingHierarchicalMultiPatchSplineSpace(
            const std::shared_ptr<const VectorConformingBasisComplex>& bc,
            const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>>& refinement_levels );
        virtual ~VectorConformingHierarchicalMultiPatchSplineSpace() = default;

        virtual const VectorConformingBasisComplex& basisComplex() const override;
        const std::shared_ptr<const VectorConformingBasisComplex>& basisComplexPtr() const { return mBasisComplex; }

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        virtual size_t numVectorComponents() const override;

        const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>>& refinementLevels() const
        {
            return mRefinementLevels;
        }

        const std::vector<std::shared_ptr<const VectorConformingHierarchicalTPSplineSpace>>& constituents() const
        {
            return mConstituents;
        }

        /// ith patch's jth local hierarchical vector function maps to this global
        /// hierarchical vector function id and orientation sign.
        const std::vector<std::vector<std::pair<FunctionId, bool>>>& functionIdMap() const { return mFuncIds; }

        private:
        const std::shared_ptr<const VectorConformingBasisComplex> mBasisComplex;
        const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>> mRefinementLevels;
        std::vector<std::shared_ptr<const VectorConformingHierarchicalTPSplineSpace>> mConstituents;
        std::vector<std::vector<std::pair<FunctionId, bool>>> mFuncIds;
        size_t mNumActiveFuncs;
    };

    VectorConformingHierarchicalMultiPatchSplineSpace
        buildHierarchicalSplineSpace(
            const std::vector<std::shared_ptr<const VectorConformingMultiPatchSplineSpace>>& refinement_levels,
            const std::vector<std::vector<topology::Cell>>& leaf_elements );

    VectorConformingHierarchicalMultiPatchSplineSpace
        buildVectorConformingHierarchicalMultiPatchSplineSpace( const HierarchicalMultiPatchSplineSpace& primal,
                                                                const ConformingType conforming_type );

    VectorConformingHierarchicalMultiPatchSplineSpace
        buildHDivHierarchicalMultiPatchSplineSpace( const HierarchicalMultiPatchSplineSpace& primal );

    VectorConformingHierarchicalMultiPatchSplineSpace
        buildHCurlHierarchicalMultiPatchSplineSpace( const HierarchicalMultiPatchSplineSpace& primal );

    HierarchicalMultiPatchSplineSpace
        buildL2HierarchicalMultiPatchSplineSpace( const VectorConformingHierarchicalMultiPatchSplineSpace& hdiv );

    std::vector<std::pair<size_t, FunctionId>>
        vectorFidsToRefinementLevelFids( const VectorConformingHierarchicalMultiPatchSplineSpace& ss );

    std::vector<VectorHierarchicalLevelFunctionId>
        patchVectorFidsToRefinementLevelFids( const VectorConformingHierarchicalMultiPatchSplineSpace& ss,
                                              const size_t patch_id );
}
