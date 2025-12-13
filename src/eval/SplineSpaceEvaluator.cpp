#include <SplineSpaceEvaluator.hpp>
#include <BasisComplex.hpp>
#include <ParentBasis.hpp>
#include <ParametricAtlas.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

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
        mParametricLengths = mSpline.basisComplex().parametricAtlas().parametricLengths( c );
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

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateParamToSpatialJacobian( const Eigen::MatrixXd& cpts ) const
    {
        const auto first_derivs = evaluateFirstDerivatives(); // Call this first to catch empty mCurrentCell
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "ParamToSpatial not supported on non-square domains" );
        return cpts( Eigen::all, mConnect ) * first_derivs * mParametricLengths.array().inverse().matrix().asDiagonal();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateParamToSpatialHessian( const Eigen::MatrixXd& cpts ) const
    {
        const auto second_derivs = evaluateSecondDerivatives(); // Call this first to catch empty mCurrentCell
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "ParamToSpatial not supported on non-square domains" );
        if( mSpline.basisComplex().parametricAtlas().cmap().dim() == 2 )
        {
            return cpts( Eigen::all, mConnect ) * second_derivs *
                   Eigen::Vector3d( mParametricLengths( 0 ) * mParametricLengths( 0 ),
                                    mParametricLengths( 0 ) * mParametricLengths( 1 ),
                                    mParametricLengths( 1 ) * mParametricLengths( 1 ) )
                       .array()
                       .inverse()
                       .matrix()
                       .asDiagonal(); // ss, ts, tt
        }
        else
        {
            return cpts( Eigen::all, mConnect ) * second_derivs *
                   Vector6d( mParametricLengths( 0 ) * mParametricLengths( 0 ),
                             mParametricLengths( 0 ) * mParametricLengths( 1 ),
                             mParametricLengths( 0 ) * mParametricLengths( 2 ),
                             mParametricLengths( 1 ) * mParametricLengths( 1 ),
                             mParametricLengths( 1 ) * mParametricLengths( 2 ),
                             mParametricLengths( 2 ) * mParametricLengths( 2 ) )
                       .array()
                       .inverse()
                       .matrix()
                       .asDiagonal(); // ss, ts, us, tt, ut, uu
        }
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateBasis() const
    {
        if( not mLocalEvals.has_value() ) throw std::runtime_error( "Must localize evaluator before evaluating" );
        return mExOp * mLocalEvals->mEvals.leftCols( mSpline.numVectorComponents() );
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateFirstDerivativesFromParamToSpatial() const
    {
        const auto first_derivs = evaluateFirstDerivatives();
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "ParamToSpatial not supported on non-square domains" );

        const size_t param_dim = mSpline.basisComplex().parametricAtlas().cmap().dim();

        const Eigen::VectorXd doubled_lengths =
            Eigen::MatrixXd( mParametricLengths.replicate( 1, param_dim ) ).transpose().reshaped();

        const Eigen::MatrixXd scaling = doubled_lengths.array().inverse().matrix().asDiagonal();

        return first_derivs * scaling;
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

    Eigen::VectorXd paramToSpatialGradDeterminant( const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd J = geom_evals.evaluateParamToSpatialJacobian( cpts );
        const Eigen::MatrixXd H = geom_evals.evaluateParamToSpatialHessian( cpts );

        if( geom_evals.splineSpace().basisComplex().parametricAtlas().cmap().dim() == 2 and cpts.rows() == 2 )
        {
            return Eigen::Vector2d( H( 0, 0 ) * J( 1, 1 ) + J( 0, 0 ) * H( 1, 1 ) - H( 1, 0 ) * J( 0, 1 ) - J( 1, 0 ) * H( 0, 1 ),
                                    H( 0, 1 ) * J( 1, 1 ) + J( 0, 0 ) * H( 1, 2 ) - H( 1, 1 ) * J( 0, 1 ) - J( 1, 0 ) * H( 0, 2 ) );
        }
        else if( geom_evals.splineSpace().basisComplex().parametricAtlas().cmap().dim() == 3 )
        {
            const double& x_s = J(0,0), x_t = J(0,1), x_u = J(0,2);
            const double& y_s = J(1,0), y_t = J(1,1), y_u = J(1,2);
            const double& z_s = J(2,0), z_t = J(2,1), z_u = J(2,2);

            const double& x_ss = H(0,0), x_st = H(0,1), x_su = H(0,2), x_tt = H(0,3), x_tu = H(0,4), x_uu = H(0,5);
            const double& y_ss = H(1,0), y_st = H(1,1), y_su = H(1,2), y_tt = H(1,3), y_tu = H(1,4), y_uu = H(1,5);
            const double& z_ss = H(2,0), z_st = H(2,1), z_su = H(2,2), z_tt = H(2,3), z_tu = H(2,4), z_uu = H(2,5);

            return Eigen::Vector3d( // Generated with mathematica
                        z_t * y_s * x_su - y_t * z_s * x_su - z_t * x_s * y_su +
                        x_t * z_s * y_su + y_t * x_s * z_su - x_t * y_s * z_su -
                        z_u * y_s * x_st + y_u * z_s * x_st + z_u * x_s * y_st -
                        x_u * z_s * y_st - y_u * x_s * z_st + x_u * y_s * z_st +
                        z_u * y_t * x_ss - y_u * z_t * x_ss - z_u * x_t * y_ss +
                        x_u * z_t * y_ss + y_u * x_t * z_ss - x_u * y_t * z_ss,

                        -z_t * y_tu * x_s + y_t * z_tu * x_s + z_u * y_tt * x_s -
                        y_u * z_tt * x_s + z_t * x_tu * y_s - x_t * z_tu * y_s -
                        z_u * x_tt * y_s + x_u * z_tt * y_s - y_t * x_tu * z_s +
                        x_t * y_tu * z_s + y_u * x_tt * z_s - x_u * y_tt * z_s +
                        z_u * y_t * x_st - y_u * z_t * x_st - z_u * x_t * y_st +
                        x_u * z_t * y_st + y_u * x_t * z_st - x_u * y_t * z_st,

                        z_uu * y_t * x_s - y_uu * z_t * x_s + z_u * y_tu * x_s -
                        y_u * z_tu * x_s - z_uu * x_t * y_s + x_uu * z_t * y_s -
                        z_u * x_tu * y_s + x_u * z_tu * y_s + y_uu * x_t * z_s -
                        x_uu * y_t * z_s + y_u * x_tu * z_s - x_u * y_tu * z_s +
                        z_u * y_t * x_su - y_u * z_t * x_su - z_u * x_t * y_su +
                        x_u * z_t * y_su + y_u * x_t * z_su - x_u * y_t * z_su
            );
        }
        else
        {
            // FIXME: Make this work on a 2d manifold in 3d
            throw std::runtime_error( "Unsupported parametric/spatial dimension for determinant gradient evaluation." );
        }
    }

    Eigen::MatrixXd piolaTransformedH1FirstDerivatives( const SplineSpaceEvaluator& scalar_evals,
                                                        const SplineSpaceEvaluator& geom_evals,
                                                        const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jac = geom_evals.evaluateJacobian( cpts );
        return scalar_evals.evaluateFirstDerivatives() *
            jac.inverse(); // transform from dxi denominator to dx
    }

    Eigen::MatrixXd piolaTransformedHCurlBasis( const SplineSpaceEvaluator& vec_evals,
                                                const SplineSpaceEvaluator& geom_evals,
                                                const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jac_inverse = geom_evals.evaluateParamToSpatialJacobian( cpts ).inverse();
        return vec_evals.evaluateBasis() * jac_inverse;
    }

    Eigen::MatrixXd doubledInverseJacobian( const Eigen::MatrixXd& jac )
    {
        const Eigen::MatrixXd inv_jac = jac.inverse();
        return kroneckerProduct( inv_jac, Eigen::MatrixXd::Identity( inv_jac.rows(), inv_jac.cols() ) );
    }

    Eigen::MatrixXd piolaTransformedHCurlFirstDerivatives( const SplineSpaceEvaluator& vec_evals,
                                                           const SplineSpaceEvaluator& geom_evals,
                                                           const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd vector_basis = vec_evals.evaluateBasis().transpose();
        const size_t param_dim = geom_evals.splineSpace().basisComplex().parametricAtlas().cmap().dim();
        const size_t spatial_dim = cpts.rows();
        const size_t n_funcs = vector_basis.cols();
        const Eigen::MatrixXd jac = geom_evals.evaluateParamToSpatialJacobian( cpts );
        const auto inverse_transpose_jacobian = jac.inverse().transpose();

        const auto modified_hessian = [&geom_evals, &cpts, &param_dim, &spatial_dim]() {
            const Eigen::MatrixXd hess = geom_evals.evaluateParamToSpatialHessian( cpts );
            Eigen::MatrixXd out( spatial_dim * param_dim, param_dim );
            if( param_dim == 2 )
                // FIXME: Make this work for 3d spatial/2d manifold
                out << hess( 0, 0 ), hess( 0, 1 ),
                       hess( 1, 0 ), hess( 1, 1 ),
                       hess( 0, 1 ), hess( 0, 2 ),
                       hess( 1, 1 ), hess( 1, 2 );
            else if( param_dim == 3 )
                out << hess( 0, 0 ), hess( 0, 1 ), hess( 0, 2 ),
                       hess( 1, 0 ), hess( 1, 1 ), hess( 1, 2 ),
                       hess( 2, 0 ), hess( 2, 1 ), hess( 2, 2 ),
                       hess( 0, 1 ), hess( 0, 3 ), hess( 0, 4 ),
                       hess( 1, 1 ), hess( 1, 3 ), hess( 1, 4 ),
                       hess( 2, 1 ), hess( 2, 3 ), hess( 2, 4 ),
                       hess( 0, 2 ), hess( 0, 4 ), hess( 0, 5 ),
                       hess( 1, 2 ), hess( 1, 4 ), hess( 1, 5 ),
                       hess( 2, 2 ), hess( 2, 4 ), hess( 2, 5 );
            else
                throw std::runtime_error( "Unsupported parametric dimension for piola transformed vector first derivatives." );
            return out;
        };

        const Eigen::MatrixXd second_term = [&]() -> Eigen::MatrixXd {
            const auto inverse_jacobian = inverse_transpose_jacobian.transpose();
            Eigen::MatrixXd multiplier = modified_hessian();

            multiplier = multiplier * inverse_jacobian;

            multiplier.topRows( spatial_dim ) = inverse_jacobian * multiplier.topRows( spatial_dim );
            multiplier.bottomRows( spatial_dim ) = inverse_jacobian * multiplier.bottomRows( spatial_dim );
            return multiplier * vector_basis;
        }();

        return ( ( inverse_transpose_jacobian *
                   vec_evals.evaluateFirstDerivativesFromParamToSpatial().transpose().reshaped( param_dim,
                                                                                                n_funcs * param_dim ) )
                     .reshaped( spatial_dim * param_dim, n_funcs ) +
                 second_term )
                   .transpose() *
               doubledInverseJacobian( jac ); // transform from ds denominator to dx
    }

    Eigen::MatrixXd piolaTransformedHDivBasis( const SplineSpaceEvaluator& vec_evals,
                                               const SplineSpaceEvaluator& geom_evals,
                                               const Eigen::MatrixXd& cpts )
    {
        const auto jac = geom_evals.evaluateParamToSpatialJacobian( cpts );
        const double det = determinant( jac );
        return 1.0 / det * vec_evals.evaluateBasis() * jac.transpose();
    }

    Eigen::MatrixXd piolaTransformedHDivFirstDerivatives( const SplineSpaceEvaluator& vec_evals,
                                                          const SplineSpaceEvaluator& geom_evals,
                                                          const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jac = geom_evals.evaluateParamToSpatialJacobian( cpts );
        const double det_inverse = 1.0 / determinant( jac );
        const Eigen::MatrixXd vector_basis = vec_evals.evaluateBasis().transpose();
        const size_t n_funcs = vector_basis.cols();
        const size_t param_dim = geom_evals.splineSpace().basisComplex().parametricAtlas().cmap().dim();
        const size_t spatial_dim = cpts.rows();

        // The product rule on the derivative of the piola transform v = ( det J )^-1 J v'
        // yields three terms from the three factors. These are the first, second, third terms here.

        const Eigen::MatrixXd first_term =
            -det_inverse * det_inverse *
            ( ( jac * vector_basis ).reshaped() * paramToSpatialGradDeterminant( geom_evals, cpts ).transpose() )
                .reshaped( spatial_dim * param_dim, n_funcs );

        const auto modified_hessian = [&geom_evals, &cpts, &param_dim, &spatial_dim]() {
            const Eigen::MatrixXd hess = geom_evals.evaluateParamToSpatialHessian( cpts );
            Eigen::MatrixXd out( spatial_dim * param_dim, param_dim );
            if( param_dim == 2 )
                // FIXME: Make this work for 3d spatial/2d manifold
                out << hess( 0, 0 ), hess( 0, 1 ),
                       hess( 1, 0 ), hess( 1, 1 ),
                       hess( 0, 1 ), hess( 0, 2 ),
                       hess( 1, 1 ), hess( 1, 2 );
            else if( param_dim == 3 )
                out << hess( 0, 0 ), hess( 0, 1 ), hess( 0, 2 ),
                       hess( 1, 0 ), hess( 1, 1 ), hess( 1, 2 ),
                       hess( 2, 0 ), hess( 2, 1 ), hess( 2, 2 ),
                       hess( 0, 1 ), hess( 0, 3 ), hess( 0, 4 ),
                       hess( 1, 1 ), hess( 1, 3 ), hess( 1, 4 ),
                       hess( 2, 1 ), hess( 2, 3 ), hess( 2, 4 ),
                       hess( 0, 2 ), hess( 0, 4 ), hess( 0, 5 ),
                       hess( 1, 2 ), hess( 1, 4 ), hess( 1, 5 ),
                       hess( 2, 2 ), hess( 2, 4 ), hess( 2, 5 );
            else
                throw std::runtime_error( "Unsupported parametric dimension for piola transformed vector first derivatives." );
            return out;
        };

        const Eigen::MatrixXd second_term = det_inverse * modified_hessian() * vector_basis;

        const Eigen::MatrixXd third_term =
            ( det_inverse * jac * vec_evals.evaluateFirstDerivativesFromParamToSpatial().transpose().reshaped( param_dim, n_funcs * param_dim ) )
                .reshaped( spatial_dim * param_dim, n_funcs );

        return ( first_term + second_term + third_term ).transpose() *
               doubledInverseJacobian( jac ); // transform from ds denominator to dx
    }

    Eigen::MatrixXd piolaTransformedL2Basis( const SplineSpaceEvaluator& bivec_evals,
                                             const SplineSpaceEvaluator& geom_evals,
                                             const Eigen::MatrixXd& cpts )
    {
        return 1.0 / determinant( geom_evals.evaluateParamToSpatialJacobian( cpts ) ) * bivec_evals.evaluateBasis();
    }
    Eigen::MatrixXd piolaTransformedL2FirstDerivatives( const SplineSpaceEvaluator& bivec_evals,
                                                        const SplineSpaceEvaluator& geom_evals,
                                                        const Eigen::MatrixXd& cpts )
    {
        const auto jac = geom_evals.evaluateParamToSpatialJacobian( cpts );
        const double det_inverse = 1.0 / determinant( jac );
        const auto jac_inverse_transpose = jac.inverse().transpose();
        return ( det_inverse * jac_inverse_transpose * bivec_evals.evaluateFirstDerivativesFromParamToSpatial().transpose() -
                 det_inverse * det_inverse * jac_inverse_transpose * paramToSpatialGradDeterminant( geom_evals, cpts ) *
                     bivec_evals.evaluateBasis().transpose() )
            .transpose() * jac.inverse(); // transform from ds denominator to dx
    }

    VertexPositionsFunc vertexPositionsFromManifold( const basis::SplineSpace& ss, const Eigen::MatrixXd& cpts )
    {
        eval::SplineSpaceEvaluator evaler( ss, 0 );
        return [evaler, cpts]( const topology::Vertex& v ) mutable -> Eigen::Vector3d {
            evaler.localizeElement( topology::Cell( v.dart(), evaler.splineSpace().basisComplex().parametricAtlas().cmap().dim() ) );
            evaler.localizePoint( evaler.splineSpace().basisComplex().parametricAtlas().parentPoint( v ) );
            return evaler.evaluateManifold( cpts );
        };
    }


    Eigen::VectorXd NURBSSpaceEvaluator::evaluateManifold( const Eigen::MatrixXd& cpts ) const
    {
        const Eigen::VectorXd bsp_eval = cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasis();
        return bsp_eval.head( bsp_eval.size() - 1 ).array() / bsp_eval( bsp_eval.size() - 1 );
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateJacobian( const Eigen::MatrixXd& cpts ) const
    {
        const Eigen::MatrixXd bsp_jac = cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateFirstDerivatives();
        const Eigen::VectorXd bsp_eval = cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasis();

        const auto DX = bsp_jac.topRows( bsp_jac.rows() - 1 );
        const auto DW = bsp_jac.bottomRows( 1 );
        const auto X = bsp_eval.head( bsp_eval.size() - 1 );
        const double w = bsp_eval( bsp_eval.size() - 1 );

        return DX / w - ( X * DW ) / ( w * w );
    }
    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateHessian( const Eigen::MatrixXd& cpts ) const
    {
        const size_t param_dim = splineSpace().basisComplex().parametricAtlas().cmap().dim();

        const Eigen::MatrixXd bsp_hess = cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateSecondDerivatives();
        const Eigen::MatrixXd bsp_jac = cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateFirstDerivatives();
        const Eigen::VectorXd bsp_eval = cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasis();

        const auto X = bsp_eval.head( bsp_eval.size() - 1 );
        const double w = bsp_eval( bsp_eval.size() - 1 );
        const auto DX = bsp_jac.topRows( bsp_jac.rows() - 1 );
        const Eigen::VectorXd DW = bsp_jac.bottomRows( 1 ).reshaped();
        const auto D2X = bsp_hess.topRows( bsp_hess.rows() - 1 );
        const Eigen::VectorXd D2W = bsp_hess.bottomRows( 1 ).reshaped();

        const double invw = 1.0 / w;
        const double invw2 = invw * invw;
        const double invw3 = invw2 * invw;

        Eigen::MatrixXd hess( D2X.rows(), D2X.cols() );

        const SmallVector<std::pair<size_t, size_t>, 6> index_map = param_dim == 2
                                                                     ? SmallVector<std::pair<size_t, size_t>, 6>{
                                                                           std::make_pair( 0, 0 ),
                                                                           std::make_pair( 0, 1 ),
                                                                           std::make_pair( 1, 1 ),
                                                                       }
                                                                     : SmallVector<std::pair<size_t, size_t>, 6>{
                                                                           std::make_pair( 0, 0 ),
                                                                           std::make_pair( 0, 1 ),
                                                                           std::make_pair( 0, 2 ),
                                                                           std::make_pair( 1, 1 ),
                                                                           std::make_pair( 1, 2 ),
                                                                           std::make_pair( 2, 2 ),
                                                                       };

        for( size_t i = 0; i < index_map.size(); ++i )
        {
            const size_t a = index_map[i].first;
            const size_t b = index_map[i].second;

            hess.col( i ) = D2X.col( i ) * invw - DX.col( a ) * DW( b ) * invw2 - DX.col( b ) * DW( a ) * invw2 -
                            X * D2W( i ) * invw2 + 2.0 * X * DW( a ) * DW( b ) * invw3;
        }

        return hess;
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateParamToSpatialJacobian( const Eigen::MatrixXd& cpts ) const
    {
        const auto jac = evaluateJacobian( cpts );
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "ParamToSpatial not supported on non-square domains" );
        return jac * mParametricLengths.array().inverse().matrix().asDiagonal();
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateParamToSpatialHessian( const Eigen::MatrixXd& cpts ) const
    {
        const auto hess = evaluateHessian( cpts );
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "ParamToSpatial not supported on non-square domains" );
        if( mSpline.basisComplex().parametricAtlas().cmap().dim() == 2 )
        {
            return hess * Eigen::Vector3d( mParametricLengths( 0 ) * mParametricLengths( 0 ),
                                           mParametricLengths( 0 ) * mParametricLengths( 1 ),
                                           mParametricLengths( 1 ) * mParametricLengths( 1 ) )
                              .array()
                              .inverse()
                              .matrix()
                              .asDiagonal(); // ss, ts, tt
        }
        else
        {
            return hess * Vector6d( mParametricLengths( 0 ) * mParametricLengths( 0 ),
                                    mParametricLengths( 0 ) * mParametricLengths( 1 ),
                                    mParametricLengths( 0 ) * mParametricLengths( 2 ),
                                    mParametricLengths( 1 ) * mParametricLengths( 1 ),
                                    mParametricLengths( 1 ) * mParametricLengths( 2 ),
                                    mParametricLengths( 2 ) * mParametricLengths( 2 ) )
                              .array()
                              .inverse()
                              .matrix()
                              .asDiagonal(); // ss, ts, us, tt, ut, uu
        }
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateBasis() const
    {
        throw std::runtime_error( "evaluateBasis not implemented for NURBSSpaceEvaluator" );
    }
    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateFirstDerivatives() const
    {
        throw std::runtime_error( "evaluateFirstDerivatives not implemented for NURBSSpaceEvaluator" );
    }
    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateSecondDerivatives() const
    {
        throw std::runtime_error( "evaluateSecondDerivatives not implemented for NURBSSpaceEvaluator" );
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateFirstDerivativesFromParamToSpatial() const
    {
        throw std::runtime_error( "evaluateFirstDerivativesFromParamToSpatial not implemented for NURBSSpaceEvaluator" );
    }
}