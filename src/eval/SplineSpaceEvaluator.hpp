#pragma once
#include <SplineSpace.hpp>
#include <ParentBasisEval.hpp>
#include <Cell.hpp>
#include <optional>
#include <CustomEigen.hpp>
#include <VertexPositionsFunc.hpp>

namespace eval
{
    class SplineSpaceEvaluator
    {
        public:
        SplineSpaceEvaluator( const basis::SplineSpace& ss, const size_t n_deriv );

        virtual ~SplineSpaceEvaluator() = default;

        void localizeElement( const topology::Cell& );
        void localizePoint( const param::ParentPoint& ppt );

        virtual Eigen::VectorXd evaluateManifold( const Eigen::MatrixXd& cpts ) const;
        virtual Eigen::MatrixXd evaluateJacobian( const Eigen::MatrixXd& cpts ) const;
        virtual Eigen::MatrixXd evaluateHessian( const Eigen::MatrixXd& cpts ) const;// FIXME: Document the returned matrix ordering. Flattening a 3-tensor into a matrix here.

        virtual Eigen::MatrixXd evaluateParamToSpatialJacobian( const Eigen::MatrixXd& cpts ) const;
        virtual Eigen::MatrixXd evaluateParamToSpatialHessian( const Eigen::MatrixXd& cpts ) const;

        virtual Eigen::MatrixXd evaluateBasis() const;
        virtual Eigen::MatrixXd evaluateFirstDerivatives() const;
        virtual Eigen::MatrixXd evaluateSecondDerivatives() const;// FIXME: Document the returned matrix ordering.

        virtual Eigen::MatrixXd evaluateFirstDerivativesFromParamToSpatial() const;

        virtual Eigen::VectorXd evaluateParametricPoint() const;

        size_t numDerivatives() const { return mNumDerivs; }

        const basis::SplineSpace& splineSpace() const { return mSpline; }

        protected:
        const size_t mNumDerivs;
        const basis::SplineSpace& mSpline;
        std::optional<topology::Cell> mCurrentCell;
        Vector6dMax mParametricLengths;
        Vector6dMax mParametricStart;
        std::vector<basis::FunctionId> mConnect;
        Eigen::MatrixXd mExOp;
        std::optional<ParentBasisEval> mLocalEvals;
        std::optional<Eigen::VectorXd> mLocalizedParentCoords;
    };

    Eigen::MatrixXd piolaTransformedH1FirstDerivatives( const SplineSpaceEvaluator& h1_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedHCurlBasis( const SplineSpaceEvaluator& vec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedHCurlFirstDerivatives( const SplineSpaceEvaluator& vec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedHDivBasis( const SplineSpaceEvaluator& vec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedHDivFirstDerivatives( const SplineSpaceEvaluator& vec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedL2Basis( const SplineSpaceEvaluator& bivec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedL2FirstDerivatives( const SplineSpaceEvaluator& bivec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );

    /// @brief Returns a VertexPositionsFunc by evaluating the manifold at each vertex.
    /// @param ss The spline space defining the manifold.
    /// @param cpts The control points defining the manifold, in a spatial_dim x n_control_pts sized matrix.
    VertexPositionsFunc vertexPositionsFromManifold( const basis::SplineSpace& ss, const Eigen::MatrixXd& cpts );


    class NURBSSpaceEvaluator : public SplineSpaceEvaluator
    {
        public:
        NURBSSpaceEvaluator( const basis::SplineSpace& ss, const size_t n_deriv )
            : SplineSpaceEvaluator( ss, n_deriv )
        {}

        virtual Eigen::VectorXd evaluateManifold( const Eigen::MatrixXd& homogeneous_cpts ) const override;
        virtual Eigen::MatrixXd evaluateJacobian( const Eigen::MatrixXd& homogeneous_cpts ) const override;
        virtual Eigen::MatrixXd evaluateHessian( const Eigen::MatrixXd& homogeneous_cpts ) const override;

        virtual Eigen::MatrixXd evaluateParamToSpatialJacobian( const Eigen::MatrixXd& homogeneous_cpts ) const override;
        virtual Eigen::MatrixXd evaluateParamToSpatialHessian( const Eigen::MatrixXd& homogeneous_cpts ) const override;

        virtual Eigen::MatrixXd evaluateBasis() const override;
        virtual Eigen::MatrixXd evaluateFirstDerivatives() const override;
        virtual Eigen::MatrixXd evaluateSecondDerivatives() const override;

        virtual Eigen::MatrixXd evaluateFirstDerivativesFromParamToSpatial() const override;
    };
}