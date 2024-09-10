#pragma once
#include <TPBasisComplex.hpp>
#include <HierarchicalTPParametricAtlas.hpp>

namespace basis
{
    /// @brief A parent basis complex over a HierarchicalTPParametricAtlas.
    ///
    class HierarchicalTPBasisComplex : public BasisComplex
    {
        public:
        HierarchicalTPBasisComplex( const std::shared_ptr<const param::HierarchicalTPParametricAtlas>& pa,
                                    const std::vector<std::shared_ptr<const TPBasisComplex>>& refinement_levels );
        virtual ~HierarchicalTPBasisComplex() = default;

        virtual const param::HierarchicalTPParametricAtlas& parametricAtlas() const override;
        const std::shared_ptr<const param::HierarchicalTPParametricAtlas>& parametricAtlasPtr() const { return mAtlas; }

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        private:
        const std::shared_ptr<const param::HierarchicalTPParametricAtlas> mAtlas;
        const std::vector<std::shared_ptr<const TPBasisComplex>> mRefinementLevels;
    };
}