#include <SplineSpaceEvaluator.hpp>
#include <BasisComplex.hpp>
#include <ParentBasis.hpp>
#include <ParametricAtlas.hpp>
#include <Eigen/Dense>
#include <string>
#include <utility>

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
        mParametricStart = mSpline.basisComplex().parametricAtlas().parametricStarts( c );
    }

    void SplineSpaceEvaluator::localizeParentPoint( const param::ParentPoint& ppt )
    {
        if( not mCurrentCell.has_value() ) throw std::runtime_error( "Must localize an element before localizing a point" );
        mLocalEvals.emplace( ParentBasisEval( mSpline.basisComplex().parentBasis( mCurrentCell.value() ), ppt, mNumDerivs ) );
        mLocalizedParentCoords.emplace( ppt.mPoint.cast<double>() );
    }

    Eigen::VectorXd SplineSpaceEvaluator::evaluateParametricPoint() const
    {
        if( not mLocalizedParentCoords.has_value() ) throw std::runtime_error( "Must localize point before evaluating parametric point" );
        const size_t param_dim = mSpline.basisComplex().parametricAtlas().cmap().dim();
        return mParametricStart.head( param_dim ) + mLocalizedParentCoords->cwiseProduct( mParametricLengths.head( param_dim ) );
    }

    Eigen::VectorXd SplineSpaceEvaluator::evaluateManifold( const Eigen::MatrixXd& cpts ) const
    {
        return cpts( Eigen::all, mConnect ) * evaluateBasisValuesAtParentPoint();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateParentToSpatialJacobian( const Eigen::MatrixXd& cpts ) const
    {
        return cpts( Eigen::all, mConnect ) * evaluateBasisFirstDerivativesWrtParentCoordinates();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateParentToSpatialHessian( const Eigen::MatrixXd& cpts ) const
    {
        return cpts( Eigen::all, mConnect ) * evaluateBasisSecondDerivativesWrtParentCoordinates();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateParametricToSpatialJacobian( const Eigen::MatrixXd& cpts ) const
    {
        const auto first_derivs = evaluateBasisFirstDerivativesWrtParentCoordinates(); // Call this first to catch empty mCurrentCell
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "Parametric-to-spatial evaluation is not supported on non-Cartesian domains" );
        return cpts( Eigen::all, mConnect ) * first_derivs * mParametricLengths.array().inverse().matrix().asDiagonal();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateParametricToSpatialHessian( const Eigen::MatrixXd& cpts ) const
    {
        const auto second_derivs = evaluateBasisSecondDerivativesWrtParentCoordinates(); // Call this first to catch empty mCurrentCell
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "Parametric-to-spatial evaluation is not supported on non-Cartesian domains" );
        if( mSpline.basisComplex().parametricAtlas().cmap().dim() == 2 )
        {
            return cpts( Eigen::all, mConnect ) * second_derivs *
                   Eigen::Vector3d( mParametricLengths( 0 ) * mParametricLengths( 0 ),
                                    mParametricLengths( 0 ) * mParametricLengths( 1 ),
                                    mParametricLengths( 1 ) * mParametricLengths( 1 ) )
                       .array()
                       .inverse()
                       .matrix()
                       .asDiagonal(); // xi-xi, xi-eta, eta-eta
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
                       .asDiagonal(); // xi-xi, xi-eta, xi-zeta, eta-eta, eta-zeta, zeta-zeta
        }
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateBasisValuesAtParentPoint() const
    {
        if( not mLocalEvals.has_value() ) throw std::runtime_error( "Must localize evaluator before evaluating" );
        return mExOp * mLocalEvals->mEvals.leftCols( mSpline.numVectorComponents() );
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateBasisFirstDerivativesWrtParametricCoordinates() const
    {
        const auto first_derivs = evaluateBasisFirstDerivativesWrtParentCoordinates();
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "Parametric-coordinate derivatives are not supported on non-Cartesian domains" );

        const size_t param_dim = mSpline.basisComplex().parametricAtlas().cmap().dim();
        const size_t n_components = mSpline.numVectorComponents();
        Eigen::VectorXd repeated_lengths( param_dim * n_components );
        for( size_t direction = 0; direction < param_dim; ++direction )
            repeated_lengths.segment( direction * n_components, n_components )
                .setConstant( mParametricLengths( direction ) );

        return first_derivs * repeated_lengths.cwiseInverse().asDiagonal();
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateBasisFirstDerivativesWrtParentCoordinates() const
    {
        if( not mLocalEvals.has_value() ) throw std::runtime_error( "Must localize evaluator before evaluating" );
        if( numDerivatives() < 1 )
            throw std::runtime_error( "Insufficient derivatives requested on SplineSpaceEvaluator construction." );
        return mExOp * mLocalEvals->mEvals.middleCols( mSpline.numVectorComponents(),
                                                      mSpline.basisComplex().parametricAtlas().cmap().dim() *
                                                          mSpline.numVectorComponents() );
    }

    Eigen::MatrixXd SplineSpaceEvaluator::evaluateBasisSecondDerivativesWrtParentCoordinates() const
    {
        if( not mLocalEvals.has_value() ) throw std::runtime_error( "Must localize evaluator before evaluating" );
        if( numDerivatives() < 2 )
            throw std::runtime_error( "Insufficient derivatives requested on SplineSpaceEvaluator construction." );
        const size_t param_dim = mSpline.basisComplex().parametricAtlas().cmap().dim();
        const size_t vec_comps = mSpline.numVectorComponents();
        return mExOp * mLocalEvals->mEvals.middleCols( vec_comps * ( 1 + param_dim ), vec_comps * param_dim * ( param_dim + 1 ) / 2 );
    }

    namespace
    {
        void requireSquareJacobian( const Eigen::MatrixXd& jac, const char* operation )
        {
            if( jac.rows() != jac.cols() )
                throw std::runtime_error( std::string( operation ) + " requires equal parametric and spatial dimensions." );
        }

        size_t symmetricDerivativeColumn( size_t first, size_t second, const size_t dim )
        {
            if( first > second ) std::swap( first, second );
            return first * dim - first * ( first + 1 ) / 2 + second;
        }

        Eigen::MatrixXd vectorBasisDerivative(
            const Eigen::MatrixXd& derivatives,
            const Eigen::Index function,
            const size_t n_components,
            const size_t parametric_dim )
        {
            Eigen::MatrixXd out( n_components, parametric_dim );
            for( size_t direction = 0; direction < parametric_dim; ++direction )
                for( size_t component = 0; component < n_components; ++component )
                    out( component, direction ) =
                        derivatives( function, direction * n_components + component );
            return out;
        }

        void storeSpatialVectorDerivative(
            Eigen::MatrixXd& output,
            const Eigen::Index function,
            const Eigen::MatrixXd& derivative )
        {
            for( Eigen::Index direction = 0; direction < derivative.cols(); ++direction )
                for( Eigen::Index component = 0; component < derivative.rows(); ++component )
                    output( function, direction * derivative.rows() + component ) =
                        derivative( component, direction );
        }
    }

    std::vector<Eigen::MatrixXd> evaluateParametricToSpatialJacobianFirstDerivatives(
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts )
    {
        const size_t parametric_dim =
            geom_evals.splineSpace().basisComplex().parametricAtlas().cmap().dim();
        const Eigen::MatrixXd hessian =
            geom_evals.evaluateParametricToSpatialHessian( cpts );
        const size_t expected_columns = parametric_dim * ( parametric_dim + 1 ) / 2;
        if( static_cast<size_t>( hessian.cols() ) != expected_columns )
            throw std::runtime_error( "Unexpected parametric-to-spatial Hessian column count." );

        std::vector<Eigen::MatrixXd> derivatives(
            parametric_dim,
            Eigen::MatrixXd( hessian.rows(), parametric_dim ) );
        for( size_t derivative_direction = 0;
             derivative_direction < parametric_dim;
             ++derivative_direction )
        {
            for( size_t jacobian_direction = 0;
                 jacobian_direction < parametric_dim;
                 ++jacobian_direction )
            {
                derivatives.at( derivative_direction ).col( jacobian_direction ) =
                    hessian.col( symmetricDerivativeColumn(
                        derivative_direction,
                        jacobian_direction,
                        parametric_dim ) );
            }
        }
        return derivatives;
    }

    Eigen::VectorXd evaluateParametricToSpatialJacobianDeterminantFirstDerivatives(
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jacobian =
            geom_evals.evaluateParametricToSpatialJacobian( cpts );
        requireSquareJacobian(
            jacobian,
            "Parametric-to-spatial Jacobian determinant differentiation" );

        const std::vector<Eigen::MatrixXd> jacobian_derivatives =
            evaluateParametricToSpatialJacobianFirstDerivatives( geom_evals, cpts );
        const Eigen::MatrixXd inverse_jacobian = jacobian.inverse();
        const double determinant = jacobian.determinant();
        Eigen::VectorXd output( jacobian_derivatives.size() );
        for( size_t direction = 0;
             direction < jacobian_derivatives.size();
             ++direction )
        {
            output( direction ) =
                determinant *
                ( inverse_jacobian * jacobian_derivatives.at( direction ) ).trace();
        }
        return output;
    }

    Eigen::MatrixXd evaluateSpatialH1BasisFirstDerivatives( const SplineSpaceEvaluator& scalar_evals,
                                                        const SplineSpaceEvaluator& geom_evals,
                                                        const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jacobian =
            geom_evals.evaluateParametricToSpatialJacobian( cpts );
        requireSquareJacobian( jacobian, "Spatial H1 derivative evaluation" );
        return scalar_evals.evaluateBasisFirstDerivativesWrtParametricCoordinates() *
               jacobian.inverse();
    }

    Eigen::MatrixXd evaluateSpatialHCurlBasisValues( const SplineSpaceEvaluator& vec_evals,
                                                const SplineSpaceEvaluator& geom_evals,
                                                const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jacobian =
            geom_evals.evaluateParametricToSpatialJacobian( cpts );
        requireSquareJacobian( jacobian, "Spatial H(curl) basis evaluation" );
        return vec_evals.evaluateBasisValuesAtParentPoint() * jacobian.inverse();
    }

    Eigen::MatrixXd evaluateSpatialHCurlBasisFirstDerivatives( const SplineSpaceEvaluator& vec_evals,
                                                           const SplineSpaceEvaluator& geom_evals,
                                                           const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jacobian =
            geom_evals.evaluateParametricToSpatialJacobian( cpts );
        requireSquareJacobian(
            jacobian,
            "Spatial H(curl) derivative evaluation" );
        const size_t parametric_dim = jacobian.cols();
        const size_t spatial_dim = jacobian.rows();
        if( vec_evals.splineSpace().numVectorComponents() != parametric_dim )
            throw std::runtime_error(
                "H(curl) component count must equal the parametric dimension." );

        const Eigen::MatrixXd inverse_jacobian = jacobian.inverse();
        const Eigen::MatrixXd inverse_transpose_jacobian =
            inverse_jacobian.transpose();
        const std::vector<Eigen::MatrixXd> jacobian_derivatives =
            evaluateParametricToSpatialJacobianFirstDerivatives(
                geom_evals,
                cpts );
        const Eigen::MatrixXd basis =
            vec_evals.evaluateBasisValuesAtParentPoint();
        const Eigen::MatrixXd basis_derivatives =
            vec_evals.evaluateBasisFirstDerivativesWrtParametricCoordinates();

        Eigen::MatrixXd output(
            basis.rows(),
            spatial_dim * spatial_dim );
        for( Eigen::Index function = 0; function < basis.rows(); ++function )
        {
            const Eigen::VectorXd parametric_value =
                basis.row( function ).transpose();
            const Eigen::MatrixXd parametric_derivative =
                vectorBasisDerivative(
                    basis_derivatives,
                    function,
                    parametric_dim,
                    parametric_dim );

            Eigen::MatrixXd spatial_derivative_wrt_parametric(
                spatial_dim,
                parametric_dim );
            for( size_t direction = 0;
                 direction < parametric_dim;
                 ++direction )
            {
                const Eigen::MatrixXd inverse_transpose_derivative =
                    -inverse_transpose_jacobian *
                    jacobian_derivatives.at( direction ).transpose() *
                    inverse_transpose_jacobian;
                spatial_derivative_wrt_parametric.col( direction ) =
                    inverse_transpose_derivative * parametric_value +
                    inverse_transpose_jacobian *
                        parametric_derivative.col( direction );
            }

            storeSpatialVectorDerivative(
                output,
                function,
                spatial_derivative_wrt_parametric * inverse_jacobian );
        }
        return output;
    }

    Eigen::MatrixXd evaluateSpatialHDivBasisValues( const SplineSpaceEvaluator& vec_evals,
                                               const SplineSpaceEvaluator& geom_evals,
                                               const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jacobian =
            geom_evals.evaluateParametricToSpatialJacobian( cpts );
        requireSquareJacobian( jacobian, "Spatial H(div) basis evaluation" );
        return vec_evals.evaluateBasisValuesAtParentPoint() *
               jacobian.transpose() / jacobian.determinant();
    }

    Eigen::MatrixXd evaluateSpatialHDivBasisFirstDerivatives( const SplineSpaceEvaluator& vec_evals,
                                                          const SplineSpaceEvaluator& geom_evals,
                                                          const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jacobian =
            geom_evals.evaluateParametricToSpatialJacobian( cpts );
        requireSquareJacobian(
            jacobian,
            "Spatial H(div) derivative evaluation" );
        const size_t parametric_dim = jacobian.cols();
        const size_t spatial_dim = jacobian.rows();
        if( vec_evals.splineSpace().numVectorComponents() != parametric_dim )
            throw std::runtime_error(
                "H(div) component count must equal the parametric dimension." );

        const double determinant = jacobian.determinant();
        const double inverse_determinant = 1.0 / determinant;
        const Eigen::MatrixXd inverse_jacobian = jacobian.inverse();
        const std::vector<Eigen::MatrixXd> jacobian_derivatives =
            evaluateParametricToSpatialJacobianFirstDerivatives(
                geom_evals,
                cpts );
        const Eigen::VectorXd determinant_derivatives =
            evaluateParametricToSpatialJacobianDeterminantFirstDerivatives(
                geom_evals,
                cpts );
        const Eigen::MatrixXd basis =
            vec_evals.evaluateBasisValuesAtParentPoint();
        const Eigen::MatrixXd basis_derivatives =
            vec_evals.evaluateBasisFirstDerivativesWrtParametricCoordinates();

        Eigen::MatrixXd output(
            basis.rows(),
            spatial_dim * spatial_dim );
        for( Eigen::Index function = 0; function < basis.rows(); ++function )
        {
            const Eigen::VectorXd parametric_value =
                basis.row( function ).transpose();
            const Eigen::MatrixXd parametric_derivative =
                vectorBasisDerivative(
                    basis_derivatives,
                    function,
                    parametric_dim,
                    parametric_dim );
            const Eigen::VectorXd jacobian_times_value =
                jacobian * parametric_value;

            Eigen::MatrixXd spatial_derivative_wrt_parametric(
                spatial_dim,
                parametric_dim );
            for( size_t direction = 0;
                 direction < parametric_dim;
                 ++direction )
            {
                spatial_derivative_wrt_parametric.col( direction ) =
                    -inverse_determinant * inverse_determinant *
                        determinant_derivatives( direction ) *
                        jacobian_times_value +
                    inverse_determinant *
                        jacobian_derivatives.at( direction ) *
                        parametric_value +
                    inverse_determinant * jacobian *
                        parametric_derivative.col( direction );
            }

            storeSpatialVectorDerivative(
                output,
                function,
                spatial_derivative_wrt_parametric * inverse_jacobian );
        }
        return output;
    }

    Eigen::MatrixXd evaluateSpatialL2BasisValues( const SplineSpaceEvaluator& l2_evals,
                                             const SplineSpaceEvaluator& geom_evals,
                                             const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jacobian =
            geom_evals.evaluateParametricToSpatialJacobian( cpts );
        requireSquareJacobian( jacobian, "Spatial L2 basis evaluation" );
        return l2_evals.evaluateBasisValuesAtParentPoint() /
               jacobian.determinant();
    }
    Eigen::MatrixXd evaluateSpatialL2BasisFirstDerivatives( const SplineSpaceEvaluator& l2_evals,
                                                        const SplineSpaceEvaluator& geom_evals,
                                                        const Eigen::MatrixXd& cpts )
    {
        const Eigen::MatrixXd jacobian =
            geom_evals.evaluateParametricToSpatialJacobian( cpts );
        requireSquareJacobian(
            jacobian,
            "Spatial L2 derivative evaluation" );
        if( l2_evals.splineSpace().numVectorComponents() != 1 )
            throw std::runtime_error(
                "L2 spatial derivative evaluation requires scalar coefficients." );

        const double inverse_determinant =
            1.0 / jacobian.determinant();
        const Eigen::VectorXd determinant_derivatives =
            evaluateParametricToSpatialJacobianDeterminantFirstDerivatives(
                geom_evals,
                cpts );
        const Eigen::MatrixXd values =
            l2_evals.evaluateBasisValuesAtParentPoint();
        const Eigen::MatrixXd parametric_derivatives =
            l2_evals.evaluateBasisFirstDerivativesWrtParametricCoordinates();
        const Eigen::MatrixXd transformed_parametric_derivatives =
            inverse_determinant * parametric_derivatives -
            inverse_determinant * inverse_determinant *
                values * determinant_derivatives.transpose();
        return transformed_parametric_derivatives * jacobian.inverse();
    }

    VertexPositionsFunc vertexPositionsFromManifold( const basis::SplineSpace& ss, const Eigen::MatrixXd& cpts )
    {
        eval::SplineSpaceEvaluator evaler( ss, 0 );
        return [evaler, cpts]( const topology::Vertex& v ) mutable -> Vector3dMax {
            evaler.localizeElement( topology::Cell( v.dart(), evaler.splineSpace().basisComplex().parametricAtlas().cmap().dim() ) );
            evaler.localizeParentPoint( evaler.splineSpace().basisComplex().parametricAtlas().parentPoint( v ) );
            return evaler.evaluateManifold( cpts );
        };
    }


    Eigen::VectorXd NURBSSpaceEvaluator::evaluateManifold( const Eigen::MatrixXd& cpts ) const
    {
        const Eigen::VectorXd bsp_eval =
            cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasisValuesAtParentPoint();
        return bsp_eval.head( bsp_eval.size() - 1 ).array() / bsp_eval( bsp_eval.size() - 1 );
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateParentToSpatialJacobian( const Eigen::MatrixXd& cpts ) const
    {
        const Eigen::MatrixXd bsp_jac =
            cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasisFirstDerivativesWrtParentCoordinates();
        const Eigen::VectorXd bsp_eval =
            cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasisValuesAtParentPoint();

        const auto DX = bsp_jac.topRows( bsp_jac.rows() - 1 );
        const auto DW = bsp_jac.bottomRows( 1 );
        const auto X = bsp_eval.head( bsp_eval.size() - 1 );
        const double w = bsp_eval( bsp_eval.size() - 1 );

        return DX / w - ( X * DW ) / ( w * w );
    }
    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateParentToSpatialHessian( const Eigen::MatrixXd& cpts ) const
    {
        const size_t param_dim = splineSpace().basisComplex().parametricAtlas().cmap().dim();

        const Eigen::MatrixXd bsp_hess =
            cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasisSecondDerivativesWrtParentCoordinates();
        const Eigen::MatrixXd bsp_jac =
            cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasisFirstDerivativesWrtParentCoordinates();
        const Eigen::VectorXd bsp_eval =
            cpts( Eigen::all, mConnect ) * SplineSpaceEvaluator::evaluateBasisValuesAtParentPoint();

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

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateParametricToSpatialJacobian( const Eigen::MatrixXd& cpts ) const
    {
        const auto jac = evaluateParentToSpatialJacobian( cpts );
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "Parametric-to-spatial evaluation is not supported on non-Cartesian domains" );
        return jac * mParametricLengths.array().inverse().matrix().asDiagonal();
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateParametricToSpatialHessian( const Eigen::MatrixXd& cpts ) const
    {
        const auto hess = evaluateParentToSpatialHessian( cpts );
        if( not param::isCartesian( mSpline.basisComplex().parametricAtlas().parentDomain( mCurrentCell.value() ) ) )
            throw std::runtime_error( "Parametric-to-spatial evaluation is not supported on non-Cartesian domains" );
        if( mSpline.basisComplex().parametricAtlas().cmap().dim() == 2 )
        {
            return hess * Eigen::Vector3d( mParametricLengths( 0 ) * mParametricLengths( 0 ),
                                           mParametricLengths( 0 ) * mParametricLengths( 1 ),
                                           mParametricLengths( 1 ) * mParametricLengths( 1 ) )
                              .array()
                              .inverse()
                              .matrix()
                              .asDiagonal(); // xi-xi, xi-eta, eta-eta
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
                              .asDiagonal(); // xi-xi, xi-eta, xi-zeta, eta-eta, eta-zeta, zeta-zeta
        }
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateBasisValuesAtParentPoint() const
    {
        throw std::runtime_error( "evaluateBasisValuesAtParentPoint not implemented for NURBSSpaceEvaluator" );
    }
    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateBasisFirstDerivativesWrtParentCoordinates() const
    {
        throw std::runtime_error( "evaluateBasisFirstDerivativesWrtParentCoordinates not implemented for NURBSSpaceEvaluator" );
    }
    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateBasisSecondDerivativesWrtParentCoordinates() const
    {
        throw std::runtime_error( "evaluateBasisSecondDerivativesWrtParentCoordinates not implemented for NURBSSpaceEvaluator" );
    }

    Eigen::MatrixXd NURBSSpaceEvaluator::evaluateBasisFirstDerivativesWrtParametricCoordinates() const
    {
        throw std::runtime_error(
            "evaluateBasisFirstDerivativesWrtParametricCoordinates not implemented for NURBSSpaceEvaluator" );
    }
}
