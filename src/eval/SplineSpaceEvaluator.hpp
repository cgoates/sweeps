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

        void localizeElement( const topology::Cell& );
        void localizePoint( const param::ParentPoint& ppt );

        Eigen::VectorXd evaluateManifold( const Eigen::MatrixXd& cpts ) const;
        Eigen::MatrixXd evaluateJacobian( const Eigen::MatrixXd& cpts ) const;
        Eigen::MatrixXd evaluateHessian( const Eigen::MatrixXd& cpts ) const;// FIXME: Document the returned matrix ordering. Flattening a 3-tensor into a matrix here.

        Eigen::MatrixXd evaluateParamToSpatialJacobian( const Eigen::MatrixXd& cpts ) const;
        Eigen::MatrixXd evaluateParamToSpatialHessian( const Eigen::MatrixXd& cpts ) const;

        Eigen::MatrixXd evaluateBasis() const;
        Eigen::MatrixXd evaluateFirstDerivatives() const;
        Eigen::MatrixXd evaluateSecondDerivatives() const;// FIXME: Document the returned matrix ordering.

        Eigen::MatrixXd evaluateFirstDerivativesFromParamToSpatial() const;

        size_t numDerivatives() const { return mNumDerivs; }

        const basis::SplineSpace& splineSpace() const { return mSpline; }

        private:
        const size_t mNumDerivs;
        const basis::SplineSpace& mSpline;
        std::optional<topology::Cell> mCurrentCell;
        Vector6dMax mParametricLengths;
        std::vector<basis::FunctionId> mConnect;
        Eigen::MatrixXd mExOp;
        std::optional<ParentBasisEval> mLocalEvals;
    };

    Eigen::MatrixXd piolaTransformedVectorBasis( const SplineSpaceEvaluator& vec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedVectorFirstDerivatives( const SplineSpaceEvaluator& vec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedBivectorBasis( const SplineSpaceEvaluator& bivec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedBivectorFirstDerivatives( const SplineSpaceEvaluator& bivec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );

    /// @brief Returns a VertexPositionsFunc by evaluating the manifold at each vertex.
    /// @param ss The spline space defining the manifold.
    /// @param cpts The control points defining the manifold, in a spatial_dim x n_control_pts sized matrix.
    VertexPositionsFunc vertexPositionsFromManifold( const basis::SplineSpace& ss, const Eigen::MatrixXd& cpts );
}