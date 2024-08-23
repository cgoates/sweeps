#include <SplineSpaceEvaluator.hpp>
#include <BasisComplex.hpp>
#include <ParentBasis.hpp>
#include <ParametricAtlas.hpp>
#include <CombinatorialMap.hpp>

namespace eval
{
    SplineSpaceEvaluator::SplineSpaceEvaluator( const basis::SplineSpace& ss, const size_t n_deriv )
        : mNumDerivs( n_deriv ), mSpline( ss )
    {}

    void SplineSpaceEvaluator::localizeElement( const topology::Cell& c )
    {
        mConnect = mSpline.connectivity( c );
        mExOp = mSpline.extractionOperator( c );
        mCurrentCell.emplace( c );
    }

    void SplineSpaceEvaluator::localizePoint( const param::ParentPoint& ppt )
    {
        if( not mCurrentCell.has_value() ) throw std::runtime_error( "Must localize cell before locaizing point" );
        mLocalEvals.emplace( ParentBasisEval( mSpline.basisComplex().parentBasis( mCurrentCell.value() ), ppt, mNumDerivs ) );
    }

    Eigen::VectorXd SplineSpaceEvaluator::evaluateManifold( const Eigen::MatrixXd& cpts ) const
    {
        return cpts( Eigen::all, mConnect ) * evaluateBasis();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateJacobian( const Eigen::MatrixXd& cpts ) const
    {
        return cpts( Eigen::all, mConnect ) * evaluateFirstDerivatives();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateHessian( const Eigen::MatrixXd& cpts ) const
    {
        return cpts( Eigen::all, mConnect ) * evaluateSecondDerivatives();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluatePiola( const Eigen::MatrixXd& cpts ) const
    {
        const auto jac = evaluateJacobian( cpts );
        const double det = jac.col( 0 ).head<3>().cross( jac.col( 1 ).head<3>() ).norm();
        return 1.0 / det * jac;
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateBasis() const
    {
        if( not mLocalEvals.has_value() ) throw std::runtime_error( "Must localize evaluator before evaluating" );
        return mExOp * mLocalEvals->mEvals.leftCols( mSpline.numVectorComponents() );
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateFirstDerivatives() const
    {
        if( not mLocalEvals.has_value() ) throw std::runtime_error( "Must localize evaluator before evaluating" );
        if( numDerivatives() < 1 )
            throw std::runtime_error( "Insufficient derivatives requested on SplineSpaceEvaluator construction." );
        return mExOp * mLocalEvals->mEvals.middleCols( mSpline.numVectorComponents(),
                                                      mSpline.basisComplex().parametricAtlas().cmap().dim() *
                                                          mSpline.numVectorComponents() );
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateSecondDerivatives() const
    {
        if( not mLocalEvals.has_value() ) throw std::runtime_error( "Must localize evaluator before evaluating" );
        if( numDerivatives() < 2 )
            throw std::runtime_error( "Insufficient derivatives requested on SplineSpaceEvaluator construction." );
        const size_t param_dim = mSpline.basisComplex().parametricAtlas().cmap().dim();
        const size_t vec_comps = mSpline.numVectorComponents();
        return mExOp * mLocalEvals->mEvals.middleCols( vec_comps * ( 1 + param_dim ), vec_comps * param_dim * ( param_dim + 1 ) / 2 );
    }

    double determinant( const Eigen::MatrixXd& jac )
    {
        if( jac.cols() == 2 and jac.rows() == 3 )
            return jac.col( 0 ).head<3>().cross( jac.col( 1 ).head<3>() ).norm();
        return jac.determinant();
    }

    Eigen::Vector2d gradDeterminant( const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts )
    {
        const Eigen::Matrix2d J = geom_evals.evaluateJacobian( cpts );
        const Eigen::MatrixXd H = geom_evals.evaluateHessian( cpts );

        return Eigen::Vector2d( H( 0, 0 ) * J( 1, 1 ) + J( 0, 0 ) * H( 1, 1 ) - H( 1, 0 ) * J( 0, 1 ) - J( 1, 0 ) * H( 0, 1 ),
                                H( 0, 1 ) * J( 1, 1 ) + J( 0, 0 ) * H( 1, 2 ) - H( 1, 1 ) * J( 0, 1 ) - J( 1, 0 ) * H( 0, 2 ) );
    }

    Eigen::MatrixXd piolaTransformedVectorBasis( const SplineSpaceEvaluator& vec_evals,
                                                 const SplineSpaceEvaluator& geom_evals,
                                                 const Eigen::MatrixXd& cpts )
    {
        const auto jac = geom_evals.evaluateJacobian( cpts );
        const double det = determinant( jac );
        return 1.0 / det * vec_evals.evaluateBasis() * jac.transpose();
    }

    Eigen::MatrixXd piolaTransformedVectorFirstDerivatives( const SplineSpaceEvaluator& vec_evals,
                                                            const SplineSpaceEvaluator& geom_evals,
                                                            const Eigen::MatrixXd& cpts )
    {
        const auto jac = geom_evals.evaluateJacobian( cpts );
        const double det_inverse = 1.0 / determinant( jac );
        const Eigen::MatrixXd vector_basis = vec_evals.evaluateBasis().transpose();
        const size_t n_funcs = vector_basis.cols();
        const size_t dim = 2;

        // The product rule on the derivative of the piola transform v = ( det J )^-1 J v'
        // yields three terms from the three factors. These are the first, second, third terms here.

        const Eigen::MatrixXd first_term =
            -det_inverse * det_inverse *
            ( ( jac * vector_basis ).reshaped() * gradDeterminant( geom_evals, cpts ).transpose() )
                .reshaped( dim * dim, n_funcs );

        const auto modified_hessian = [&geom_evals, &cpts, &dim]() {
            const Eigen::MatrixXd hess = geom_evals.evaluateHessian( cpts );
            Eigen::MatrixXd out( dim * dim, dim );
            out << hess( 0, 0 ), hess( 0, 1 ), hess( 1, 0 ), hess( 1, 1 ), hess( 0, 1 ), hess( 0, 2 ), hess( 1, 1 ),
                hess( 1, 2 );
            return out;
        };

        const Eigen::MatrixXd second_term = det_inverse * modified_hessian() * vector_basis;

        const Eigen::MatrixXd third_term =
            ( det_inverse * jac * vec_evals.evaluateFirstDerivatives().transpose().reshaped( dim, n_funcs * dim ) )
                .reshaped( dim * dim, n_funcs );

        return ( first_term + second_term + third_term ).transpose();
    }

    Eigen::VectorXd piolaTransformedBivectorBasis( const SplineSpaceEvaluator& bivec_evals,
                                                   const SplineSpaceEvaluator& geom_evals,
                                                   const Eigen::MatrixXd& cpts )
    {
        return 1.0 / determinant( geom_evals.evaluateJacobian( cpts ) ) * bivec_evals.evaluateBasis();
    }
    Eigen::MatrixXd piolaTransformedBivectorFirstDerivatives( const SplineSpaceEvaluator& bivec_evals,
                                                              const SplineSpaceEvaluator& geom_evals,
                                                              const Eigen::MatrixXd& cpts )
    {
        const auto jac = geom_evals.evaluateJacobian( cpts );
        const double det_inverse = 1.0 / determinant( jac );
        const auto jac_inverse_transpose = jac.inverse().transpose();
        return ( det_inverse * jac_inverse_transpose * bivec_evals.evaluateFirstDerivatives().transpose() -
                 det_inverse * det_inverse * jac_inverse_transpose * gradDeterminant( geom_evals, cpts ) *
                     bivec_evals.evaluateBasis().transpose() )
            .transpose();
    }
}