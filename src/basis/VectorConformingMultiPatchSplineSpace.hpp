#pragma once
#include <VectorConformingTPSplineSpace.hpp>
#include <VectorConformingBasisComplex.hpp>
#include <MultiPatchParametricAtlas.hpp>

namespace basis
{
    class VectorConformingMultiPatchSplineSpace : public SplineSpace
    {
        public:
        VectorConformingMultiPatchSplineSpace( const std::shared_ptr<const VectorConformingBasisComplex>& bc,
                                               const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>& constituents );

        virtual ~VectorConformingMultiPatchSplineSpace() = default;

        virtual const VectorConformingBasisComplex& basisComplex() const override;
        const std::shared_ptr<const VectorConformingBasisComplex>& basisComplexPtr() const { return mBasisComplex; }

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        virtual size_t numVectorComponents() const override;

        const std::vector<std::vector<std::pair<FunctionId, bool>>>& functionIdMap() const { return mFuncIds; }

        const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>& subSpaces() const { return mSubSpaces; }

        private:
        const std::shared_ptr<const VectorConformingBasisComplex> mBasisComplex;
        const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>> mSubSpaces;
        const param::MultiPatchParametricAtlas& mParametricAtlas;//To get to the underlying multipatch structure
        std::vector<std::vector<std::pair<FunctionId, bool>>> mFuncIds;
        size_t mNumFunctions;
    };
} // namespace basis