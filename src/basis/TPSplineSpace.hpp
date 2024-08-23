#pragma once
#include <BSplineSpace1d.hpp>
#include <TPBasisComplex.hpp>

namespace basis
{
    class TPSplineSpace : public SplineSpace
    {
        public:
        TPSplineSpace( const std::shared_ptr<const TPBasisComplex>& bc,
                       const std::shared_ptr<const SplineSpace>& source,
                       const std::shared_ptr<const BSplineSpace1d>& line );
        virtual ~TPSplineSpace() = default;

        virtual const TPBasisComplex& basisComplex() const override;
        const std::shared_ptr<const TPBasisComplex>& basisComplexPtr() const { return mBasisComplex; }

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        const BSplineSpace1d& line() const { return *mLine; }
        const SplineSpace& source() const { return *mSource; }

        const std::shared_ptr<const BSplineSpace1d>& linePtr() const { return mLine; }
        const std::shared_ptr<const SplineSpace>& sourcePtr() const { return mSource; }

        private:

        FunctionId flatten( const FunctionId& source_fid, const FunctionId& line_fid ) const;

        const std::shared_ptr<const TPBasisComplex> mBasisComplex;
        const std::shared_ptr<const SplineSpace> mSource;
        const std::shared_ptr<const BSplineSpace1d> mLine;
    };

    SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> tensorProductComponentSplines( const TPSplineSpace& ss );

    TPSplineSpace buildBSpline( const SmallVector<KnotVector, 3>& kvs, const SmallVector<size_t, 3>& degrees );
}