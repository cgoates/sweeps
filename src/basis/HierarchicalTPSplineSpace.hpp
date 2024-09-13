#pragma once
#include <TPSplineSpace.hpp>
#include <HierarchicalTPBasisComplex.hpp>

namespace basis
{
    class HierarchicalTPSplineSpace : public SplineSpace
    {
        public:
        HierarchicalTPSplineSpace( const std::shared_ptr<const HierarchicalTPBasisComplex>& bc,
                                   const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels,
                                   const std::vector<std::vector<FunctionId>>& active_funcs );
        virtual ~HierarchicalTPSplineSpace() = default;

        virtual const HierarchicalTPBasisComplex& basisComplex() const override;
        const std::shared_ptr<const HierarchicalTPBasisComplex>& basisComplexPtr() const;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        private:

        const std::shared_ptr<const HierarchicalTPBasisComplex> mBasisComplex;
        const std::vector<std::shared_ptr<const TPSplineSpace>> mRefinementLevels;
        std::vector<Eigen::SparseMatrix<double>> mLevelExtractionOps;
    };
}