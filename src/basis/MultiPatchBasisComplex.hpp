#pragma once
#include <TPBasisComplex.hpp>
#include <MultiPatchParametricAtlas.hpp>

namespace basis
{
    class MultiPatchBasisComplex : public BasisComplex
    {
        public:
        MultiPatchBasisComplex( const std::shared_ptr<const param::MultiPatchParametricAtlas>& p_atlas,
                                   const std::vector<std::shared_ptr<const TPBasisComplex>>& constituents );

        virtual ~MultiPatchBasisComplex() = default;

        virtual const param::MultiPatchParametricAtlas& parametricAtlas() const override;

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        const std::vector<std::shared_ptr<const TPBasisComplex>>& constituents() const { return mSubComplexes; }

        private:
        const std::shared_ptr<const param::MultiPatchParametricAtlas> mAtlas;
        const std::vector<std::shared_ptr<const TPBasisComplex>> mSubComplexes;
    };
}