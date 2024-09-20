#pragma once
#include <TPSplineSpace.hpp>
#include <HierarchicalTPBasisComplex.hpp>

namespace basis
{
    class HierarchicalTPSplineSpace : public SplineSpace
    {
        public:
        HierarchicalTPSplineSpace( const std::shared_ptr<const HierarchicalTPBasisComplex>& bc,
                                   const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels );
        virtual ~HierarchicalTPSplineSpace() = default;

        virtual const HierarchicalTPBasisComplex& basisComplex() const override;
        const std::shared_ptr<const HierarchicalTPBasisComplex>& basisComplexPtr() const;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        const std::vector<std::shared_ptr<const TPSplineSpace>>& refinementLevels() const { return mRefinementLevels; }

        const std::vector<std::vector<FunctionId>>& activeFunctions() const { return mActiveFunctions; }

        const std::vector<Eigen::SparseMatrix<double>>& levelExtractionOperators() const { return mLevelExtractionOps; }

        private:

        const std::shared_ptr<const HierarchicalTPBasisComplex> mBasisComplex;
        const std::vector<std::shared_ptr<const TPSplineSpace>> mRefinementLevels;
        std::vector<Eigen::SparseMatrix<double>> mLevelExtractionOps;
        std::vector<std::vector<FunctionId>> mActiveFunctions;
    };

    HierarchicalTPSplineSpace
        buildHierarchicalSplineSpace( const std::vector<std::shared_ptr<const TPSplineSpace>>& refinement_levels,
                                      const std::vector<std::vector<topology::Cell>>& leaf_elements );

    /// @brief Builds the operator that transforms coefficients from the unrefined space to the hierarchical space.
    /// @param hier_ss The hierarchical space to transform the coefficients to.  Coefficients should correspond to
    ///                the spline space hier_ss.refinementLevels().at( 0 ).
    /// @return An operator to transform the coefficients.
    Eigen::SparseMatrix<double> prolongationOperator( const HierarchicalTPSplineSpace& hier_ss );
}