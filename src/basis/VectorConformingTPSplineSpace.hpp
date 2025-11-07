#pragma once
#include <SplineSpace.hpp>
#include <VectorConformingBasisComplex.hpp>
#include <TPSplineSpace.hpp>

namespace basis
{
    class VectorConformingTPSplineSpace : public SplineSpace
    {
        public:
        VectorConformingTPSplineSpace( const std::shared_ptr<const VectorConformingBasisComplex>& bc,
                                       const TPSplineSpace& primal_basis );
        virtual ~VectorConformingTPSplineSpace() = default;

        virtual const VectorConformingBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        virtual size_t numVectorComponents() const override;

        const SmallVector<std::shared_ptr<const BSplineSpace1d>, 3>& reducedDegree1dBases() const
        {
            return mReducedDegree1dBases;
        }

        const SmallVector<std::shared_ptr<const TPSplineSpace>, 3>& scalarTPBases() const { return mScalarTPBases; }

        const ConformingType& conformingType() const { return mConformingType; }

        private:
        const std::shared_ptr<const VectorConformingBasisComplex> mBasisComplex;
        SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> mReducedDegree1dBases;
        SmallVector<std::shared_ptr<const TPSplineSpace>, 3> mScalarTPBases;
        const ConformingType mConformingType;
    };
}