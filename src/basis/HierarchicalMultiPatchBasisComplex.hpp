#pragma once
#include <MultiPatchBasisComplex.hpp>
#include <HierarchicalMultiPatchParametricAtlas.hpp>
#include <HierarchicalTPBasisComplex.hpp>

namespace basis
{
    /// @brief A parent basis complex over a HierarchicalTPParametricAtlas.
    ///
    class HierarchicalMultiPatchBasisComplex : public BasisComplex
    {
        public:
        HierarchicalMultiPatchBasisComplex( const std::shared_ptr<const param::HierarchicalMultiPatchParametricAtlas>& pa,
                                    const std::vector<std::shared_ptr<const MultiPatchBasisComplex>>& refinement_levels );
        virtual ~HierarchicalMultiPatchBasisComplex() = default;

        virtual const param::HierarchicalMultiPatchParametricAtlas& parametricAtlas() const override { return *mAtlas; }
        const std::shared_ptr<const param::HierarchicalMultiPatchParametricAtlas>& parametricAtlasPtr() const { return mAtlas; }

        const std::vector<std::shared_ptr<const MultiPatchBasisComplex>>& refinementLevels() const { return mRefinementLevels; }
        const std::vector<std::shared_ptr<const HierarchicalTPBasisComplex>>& constituents() const { return mConstituents; }

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        private:
        const std::shared_ptr<const param::HierarchicalMultiPatchParametricAtlas> mAtlas;
        const std::vector<std::shared_ptr<const MultiPatchBasisComplex>> mRefinementLevels;
        const std::vector<std::shared_ptr<const HierarchicalTPBasisComplex>> mConstituents;
    };
}