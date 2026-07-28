#include <catch2/catch_test_macros.hpp>
#include <BSplineSpace1d.hpp>
#include <BasisComplex1d.hpp>
#include <ParametricAtlas1d.hpp>
#include <CombinatorialMap1d.hpp>
#include <KnotVector.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <TPSplineSpace.hpp>
#include <VTKOutput.hpp>
#include <IndexOperations.hpp>
#include <VectorConformingBasisComplex.hpp>
#include <VectorConformingTPSplineSpace.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <VTKOutput.hpp>
#include <Eigen/Dense>
#include <array>

using namespace basis;
using namespace param;
using namespace topology;

constexpr bool LOG_EVALS = false;

TEST_CASE( "TP Spline space evaluation" )
{
    // Input info
    const KnotVector kv1( {0,0,0,0,1,2,3,4,4,4,4}, 1e-10 );
    const KnotVector kv2( {0,0,0,1,2,2,2}, 1e-10 );
    const size_t degree1 = 3;
    const size_t degree2 = 2;

    // Built from input info
    const auto cmap1 = std::make_shared<const CombinatorialMap1d>( numElements( kv1 ) );
    const auto cmap2 = std::make_shared<const CombinatorialMap1d>( numElements( kv2 ) );
    const auto param1 = std::make_shared<const ParametricAtlas1d>( cmap1, parametricLengths( kv1 ) );
    const auto param2 = std::make_shared<const ParametricAtlas1d>( cmap2, parametricLengths( kv2 ) );
    const auto bc1 = std::make_shared<const BasisComplex1d>( param1, degree1 );
    const auto bc2 = std::make_shared<const BasisComplex1d>( param2, degree2 );
    const auto ss1 = std::make_shared<const BSplineSpace1d>( bc1, kv1 );
    const auto ss2 = std::make_shared<const BSplineSpace1d>( bc2, kv2 );

    const auto cmap_2d = std::make_shared<const TPCombinatorialMap>( cmap1, cmap2 );
    const auto param_2d = std::make_shared<const TPParametricAtlas>( cmap_2d, param1, param2 );
    const auto bc_2d = std::make_shared<const TPBasisComplex>( param_2d, bc1, bc2 );
    const TPSplineSpace ss_2d( bc_2d, ss1, ss2 );

    eval::SplineSpaceEvaluator primal_evals( ss_2d, 1 );

    iterateCellsWhile( *cmap_2d, 1, [&]( const topology::Edge& e ) {
        const auto maybe_opp_e = phi( *cmap_2d, 2, e.dart() );
        if( maybe_opp_e.has_value() )
        {
            const topology::Face f( e.dart() );
            const topology::Face f_opp( maybe_opp_e.value() );
            const ParentPoint ppt =
                pointOnBoundary( param_2d->parentDomain( f ), parentDomainBoundary( *param_2d, e ) );
            primal_evals.localizeElement( f );
            primal_evals.localizeParentPoint( ppt );

            const auto conn = ss_2d.connectivity( f );
            const auto basis = primal_evals.evaluateBasisValuesAtParentPoint();
            const auto derv = primal_evals.evaluateBasisFirstDerivativesWrtParentCoordinates();

            const ParentPoint ppt_opp =
                pointOnBoundary( param_2d->parentDomain( f_opp ),
                                 parentDomainBoundary( *param_2d, topology::Edge( maybe_opp_e.value() ) ) );
            primal_evals.localizeElement( f_opp );
            primal_evals.localizeParentPoint( ppt_opp );

            const auto conn_opp = ss_2d.connectivity( f_opp );
            const auto basis_opp = primal_evals.evaluateBasisValuesAtParentPoint();
            const auto derv_opp = primal_evals.evaluateBasisFirstDerivativesWrtParentCoordinates();

            for( size_t i = 0; i < conn.size(); i++ )
            {
                const auto it = std::find( conn_opp.begin(), conn_opp.end(), conn.at( i ) );
                if( it != conn_opp.end() )
                {
                    const size_t j = std::distance( conn_opp.begin(), it );

                    CHECK( util::equals( basis.row( i ), basis_opp.row( j ), 1e-9 ) );
                    CHECK( util::equals( derv.row( i ), derv_opp.row( j ), 1e-9 ) );

                    if( not util::equals( derv.row( i ), derv_opp.row( j ), 1e-9 ) )
                    {
                        std::cout << e.dart() << ", " << conn.at( i ) << ":\n";
                        std::cout << basis.row( i ) << " vs " << basis_opp.row( j ) << std::endl;
                        std::cout << derv.row( i ) << " vs " << derv_opp.row( j ) << std::endl;
                        std::cout << "----------\n";
                    }
                }
            }
        }
        return true;
    } );
}

TEST_CASE( "Parent and parametric evaluator coordinate contracts" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s( { 2, 2, 5, 5 }, ptol );
    const KnotVector kv_t( { 7, 7, 9, 9 }, ptol );
    const TPSplineSpace h1 = buildBSpline( { kv_s, kv_t }, { 1, 1 } );

    const topology::Face elem( topology::Dart( 0 ) );
    const Eigen::Vector2d parent_coords( 0.37, 0.61 );
    const param::ParentPoint parent_point(
        param::cubeDomain( 2 ), parent_coords, { false, false, false, false } );

    eval::SplineSpaceEvaluator h1_evals( h1, 2 );
    h1_evals.localizeElement( elem );
    h1_evals.localizeParentPoint( parent_point );

    const auto& atlas = h1.basisComplex().parametricAtlas();
    const Eigen::Vector2d parametric_start = atlas.parametricStarts( elem ).head<2>();
    const Eigen::Vector2d parametric_lengths = atlas.parametricLengths( elem ).head<2>();
    CHECK( ( parametric_lengths - Eigen::Vector2d( 3.0, 2.0 ) ).norm() < 1e-12 );
    const Eigen::Vector2d expected_parametric_point =
        parametric_start + parent_coords.cwiseProduct( parametric_lengths );
    CHECK( ( h1_evals.evaluateParametricPoint() - expected_parametric_point ).norm() < 1e-12 );

    Eigen::MatrixXd geometry_control_points( 2, h1.numFunctions() );
    for( Eigen::Index i = 0; i < geometry_control_points.cols(); ++i )
    {
        geometry_control_points( 0, i ) = static_cast<double>( i ) + 0.25 * static_cast<double>( i * i );
        geometry_control_points( 1, i ) = 0.5 * static_cast<double>( i ) + 0.1 * static_cast<double>( i * i );
    }

    const Eigen::MatrixXd parent_to_spatial =
        h1_evals.evaluateParentToSpatialJacobian( geometry_control_points );
    const Eigen::MatrixXd parametric_to_spatial =
        h1_evals.evaluateParametricToSpatialJacobian( geometry_control_points );
    CHECK( ( parent_to_spatial -
             parametric_to_spatial * parametric_lengths.asDiagonal() )
               .norm() < 1e-12 );

    const Eigen::MatrixXd parent_to_spatial_hessian =
        h1_evals.evaluateParentToSpatialHessian( geometry_control_points );
    const Eigen::MatrixXd parametric_to_spatial_hessian =
        h1_evals.evaluateParametricToSpatialHessian( geometry_control_points );
    const Eigen::Vector3d hessian_scaling(
        parametric_lengths( 0 ) * parametric_lengths( 0 ),
        parametric_lengths( 0 ) * parametric_lengths( 1 ),
        parametric_lengths( 1 ) * parametric_lengths( 1 ) );
    CHECK( ( parent_to_spatial_hessian -
             parametric_to_spatial_hessian * hessian_scaling.asDiagonal() )
               .norm() < 1e-12 );

    const Eigen::MatrixXd parent_scalar_derivatives =
        h1_evals.evaluateBasisFirstDerivativesWrtParentCoordinates();
    const Eigen::MatrixXd parametric_scalar_derivatives =
        h1_evals.evaluateBasisFirstDerivativesWrtParametricCoordinates();
    CHECK( ( parametric_scalar_derivatives -
             parent_scalar_derivatives * parametric_lengths.cwiseInverse().asDiagonal() )
               .norm() < 1e-12 );

    const auto hdiv_basis_complex =
        std::make_shared<const VectorConformingBasisComplex>( h1.basisComplexPtr() );
    const VectorConformingTPSplineSpace hdiv( hdiv_basis_complex, h1 );
    eval::SplineSpaceEvaluator hdiv_evals( hdiv, 1 );
    hdiv_evals.localizeElement( elem );
    hdiv_evals.localizeParentPoint( parent_point );

    const Eigen::MatrixXd parent_vector_derivatives =
        hdiv_evals.evaluateBasisFirstDerivativesWrtParentCoordinates();
    const Eigen::MatrixXd parametric_vector_derivatives =
        hdiv_evals.evaluateBasisFirstDerivativesWrtParametricCoordinates();
    const Eigen::Vector4d repeated_inverse_lengths(
        1.0 / parametric_lengths( 0 ),
        1.0 / parametric_lengths( 0 ),
        1.0 / parametric_lengths( 1 ),
        1.0 / parametric_lengths( 1 ) );
    CHECK( ( parametric_vector_derivatives -
             parent_vector_derivatives * repeated_inverse_lengths.asDiagonal() )
               .norm() < 1e-12 );
}

namespace
{
    param::ParentPoint parentPoint2d( const Eigen::Vector2d& coordinates )
    {
        return param::ParentPoint(
            param::cubeDomain( 2 ),
            coordinates,
            { false, false, false, false } );
    }

    void localizeAt(
        eval::SplineSpaceEvaluator& basis_evaluator,
        eval::SplineSpaceEvaluator& geometry_evaluator,
        const Eigen::Vector2d& coordinates )
    {
        const param::ParentPoint point = parentPoint2d( coordinates );
        basis_evaluator.localizeParentPoint( point );
        geometry_evaluator.localizeParentPoint( point );
    }

    struct SpatialTestSpaces
    {
        SpatialTestSpaces()
            : h1( buildBSpline(
                  { KnotVector( { 0, 0, 0, 3, 3, 3 }, 1e-10 ),
                    KnotVector( { 0, 0, 0, 2, 2, 2 }, 1e-10 ) },
                  { 2, 2 } ) ),
              hcurl_complex( std::make_shared<const VectorConformingBasisComplex>(
                  h1.basisComplexPtr(),
                  ConformingType::Curl ) ),
              hcurl( hcurl_complex, h1 ),
              hdiv_complex( std::make_shared<const VectorConformingBasisComplex>(
                  h1.basisComplexPtr(),
                  ConformingType::Divergence ) ),
              hdiv( hdiv_complex, h1 ),
              l2_complex( std::make_shared<const TPBasisComplex>(
                  h1.basisComplex().parametricAtlasPtr(),
                  hdiv.reducedDegree1dBases().at( 0 )->basisComplexPtr(),
                  hdiv.reducedDegree1dBases().at( 1 )->basisComplexPtr() ) ),
              l2( l2_complex,
                  hdiv.reducedDegree1dBases().at( 0 ),
                  hdiv.reducedDegree1dBases().at( 1 ) ),
              geometry_control_points( 2, h1.numFunctions() )
        {
            for( size_t j = 0; j < 3; ++j )
            {
                for( size_t i = 0; i < 3; ++i )
                {
                    const double s_control = 0.5 * static_cast<double>( i );
                    const double t_control = 0.5 * static_cast<double>( j );
                    const Eigen::Index column = static_cast<Eigen::Index>( 3 * j + i );
                    geometry_control_points( 0, column ) =
                        s_control + 0.2 * s_control * t_control;
                    geometry_control_points( 1, column ) =
                        t_control + ( i == 2 ? 0.1 : 0.0 );
                }
            }
        }

        TPSplineSpace h1;
        std::shared_ptr<const VectorConformingBasisComplex> hcurl_complex;
        VectorConformingTPSplineSpace hcurl;
        std::shared_ptr<const VectorConformingBasisComplex> hdiv_complex;
        VectorConformingTPSplineSpace hdiv;
        std::shared_ptr<const TPBasisComplex> l2_complex;
        TPSplineSpace l2;
        Eigen::MatrixXd geometry_control_points;
        topology::Face element{ topology::Dart( 0 ) };
    };

    template <typename ValueEvaluator>
    Eigen::MatrixXd finiteDifferenceSpatialDerivatives(
        eval::SplineSpaceEvaluator& basis_evaluator,
        eval::SplineSpaceEvaluator& geometry_evaluator,
        const Eigen::MatrixXd& geometry_control_points,
        const Eigen::Vector2d& parent_coordinates,
        const Eigen::Vector2d& parametric_lengths,
        ValueEvaluator&& evaluate_values )
    {
        constexpr double epsilon = 1e-6;
        std::array<Eigen::MatrixXd, 2> parametric_derivatives;
        for( size_t direction = 0; direction < 2; ++direction )
        {
            Eigen::Vector2d parent_step = Eigen::Vector2d::Zero();
            parent_step( direction ) = epsilon / parametric_lengths( direction );
            localizeAt(
                basis_evaluator,
                geometry_evaluator,
                parent_coordinates + parent_step );
            const Eigen::MatrixXd plus = evaluate_values();
            localizeAt(
                basis_evaluator,
                geometry_evaluator,
                parent_coordinates - parent_step );
            const Eigen::MatrixXd minus = evaluate_values();
            parametric_derivatives.at( direction ) =
                ( plus - minus ) / ( 2.0 * epsilon );
        }

        localizeAt(
            basis_evaluator,
            geometry_evaluator,
            parent_coordinates );
        const Eigen::MatrixXd inverse_jacobian =
            geometry_evaluator
                .evaluateParametricToSpatialJacobian(
                    geometry_control_points )
                .inverse();
        const Eigen::Index n_functions =
            parametric_derivatives.front().rows();
        const Eigen::Index n_components =
            parametric_derivatives.front().cols();
        Eigen::MatrixXd output(
            n_functions,
            n_components * inverse_jacobian.cols() );
        for( Eigen::Index function = 0; function < n_functions; ++function )
        {
            Eigen::MatrixXd derivative_wrt_parametric( n_components, 2 );
            for( size_t direction = 0; direction < 2; ++direction )
                derivative_wrt_parametric.col( direction ) =
                    parametric_derivatives.at( direction )
                        .row( function )
                        .transpose();
            const Eigen::MatrixXd derivative_wrt_spatial =
                derivative_wrt_parametric * inverse_jacobian;
            for( Eigen::Index direction = 0;
                 direction < derivative_wrt_spatial.cols();
                 ++direction )
            {
                for( Eigen::Index component = 0;
                     component < derivative_wrt_spatial.rows();
                     ++component )
                {
                    output(
                        function,
                        direction * n_components + component ) =
                        derivative_wrt_spatial( component, direction );
                }
            }
        }
        return output;
    }
}

TEST_CASE( "Spatial transforms agree through parent and parametric maps" )
{
    SpatialTestSpaces spaces;
    eval::SplineSpaceEvaluator geometry_evaluator( spaces.h1, 2 );
    eval::SplineSpaceEvaluator h1_evaluator( spaces.h1, 1 );
    eval::SplineSpaceEvaluator hcurl_evaluator( spaces.hcurl, 1 );
    eval::SplineSpaceEvaluator hdiv_evaluator( spaces.hdiv, 1 );
    eval::SplineSpaceEvaluator l2_evaluator( spaces.l2, 1 );
    for( eval::SplineSpaceEvaluator* evaluator :
         { &geometry_evaluator,
           &h1_evaluator,
           &hcurl_evaluator,
           &hdiv_evaluator,
           &l2_evaluator } )
        evaluator->localizeElement( spaces.element );

    const Eigen::Vector2d parent_coordinates( 0.37, 0.41 );
    localizeAt(
        h1_evaluator,
        geometry_evaluator,
        parent_coordinates );
    hcurl_evaluator.localizeParentPoint( parentPoint2d( parent_coordinates ) );
    hdiv_evaluator.localizeParentPoint( parentPoint2d( parent_coordinates ) );
    l2_evaluator.localizeParentPoint( parentPoint2d( parent_coordinates ) );

    const Eigen::Vector2d lengths =
        spaces.h1.basisComplex()
            .parametricAtlas()
            .parametricLengths( spaces.element )
            .head<2>();
    const Eigen::Matrix2d parent_to_parametric =
        lengths.asDiagonal();
    const double parent_to_parametric_determinant =
        parent_to_parametric.determinant();
    const Eigen::MatrixXd parent_jacobian =
        geometry_evaluator.evaluateParentToSpatialJacobian(
            spaces.geometry_control_points );

    const Eigen::MatrixXd h1_via_parent =
        h1_evaluator.evaluateBasisFirstDerivativesWrtParentCoordinates() *
        parent_jacobian.inverse();
    CHECK(
        ( h1_via_parent -
          eval::evaluateSpatialH1BasisFirstDerivatives(
              h1_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) )
                .norm() < 1e-12 );

    const Eigen::MatrixXd hcurl_parent_components =
        hcurl_evaluator.evaluateBasisValuesAtParentPoint() *
        parent_to_parametric;
    const Eigen::MatrixXd hcurl_via_parent =
        hcurl_parent_components * parent_jacobian.inverse();
    CHECK(
        ( hcurl_via_parent -
          eval::evaluateSpatialHCurlBasisValues(
              hcurl_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) )
                .norm() < 1e-12 );

    const Eigen::MatrixXd hdiv_parent_components =
        parent_to_parametric_determinant *
        hdiv_evaluator.evaluateBasisValuesAtParentPoint() *
        parent_to_parametric.inverse();
    const Eigen::MatrixXd hdiv_via_parent =
        hdiv_parent_components * parent_jacobian.transpose() /
        parent_jacobian.determinant();
    CHECK(
        ( hdiv_via_parent -
          eval::evaluateSpatialHDivBasisValues(
              hdiv_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) )
                .norm() < 1e-12 );

    const Eigen::MatrixXd l2_via_parent =
        parent_to_parametric_determinant *
        l2_evaluator.evaluateBasisValuesAtParentPoint() /
        parent_jacobian.determinant();
    CHECK(
        ( l2_via_parent -
          eval::evaluateSpatialL2BasisValues(
              l2_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) )
                .norm() < 1e-12 );
}

TEST_CASE( "Parametric geometry derivatives agree with finite differences" )
{
    SpatialTestSpaces spaces;
    eval::SplineSpaceEvaluator geometry_evaluator( spaces.h1, 2 );
    geometry_evaluator.localizeElement( spaces.element );
    const Eigen::Vector2d parent_coordinates( 0.37, 0.41 );
    const Eigen::Vector2d lengths =
        spaces.h1.basisComplex()
            .parametricAtlas()
            .parametricLengths( spaces.element )
            .head<2>();
    geometry_evaluator.localizeParentPoint(
        parentPoint2d( parent_coordinates ) );

    const std::vector<Eigen::MatrixXd> analytic_jacobian_derivatives =
        eval::evaluateParametricToSpatialJacobianFirstDerivatives(
            geometry_evaluator,
            spaces.geometry_control_points );
    const Eigen::VectorXd analytic_determinant_derivatives =
        eval::evaluateParametricToSpatialJacobianDeterminantFirstDerivatives(
            geometry_evaluator,
            spaces.geometry_control_points );

    constexpr double epsilon = 1e-6;
    for( size_t direction = 0; direction < 2; ++direction )
    {
        Eigen::Vector2d parent_step = Eigen::Vector2d::Zero();
        parent_step( direction ) = epsilon / lengths( direction );
        geometry_evaluator.localizeParentPoint(
            parentPoint2d( parent_coordinates + parent_step ) );
        const Eigen::MatrixXd plus_jacobian =
            geometry_evaluator.evaluateParametricToSpatialJacobian(
                spaces.geometry_control_points );
        geometry_evaluator.localizeParentPoint(
            parentPoint2d( parent_coordinates - parent_step ) );
        const Eigen::MatrixXd minus_jacobian =
            geometry_evaluator.evaluateParametricToSpatialJacobian(
                spaces.geometry_control_points );

        const Eigen::MatrixXd finite_difference_jacobian =
            ( plus_jacobian - minus_jacobian ) / ( 2.0 * epsilon );
        const double finite_difference_determinant =
            ( plus_jacobian.determinant() -
              minus_jacobian.determinant() ) /
            ( 2.0 * epsilon );
        CHECK(
            ( analytic_jacobian_derivatives.at( direction ) -
              finite_difference_jacobian )
                    .norm() < 1e-8 );
        CHECK(
            std::abs(
                analytic_determinant_derivatives( direction ) -
                finite_difference_determinant ) < 1e-8 );
    }
}

TEST_CASE( "Spatial basis derivatives agree with finite differences" )
{
    SpatialTestSpaces spaces;
    eval::SplineSpaceEvaluator geometry_evaluator( spaces.h1, 2 );
    eval::SplineSpaceEvaluator h1_evaluator( spaces.h1, 1 );
    eval::SplineSpaceEvaluator hcurl_evaluator( spaces.hcurl, 1 );
    eval::SplineSpaceEvaluator hdiv_evaluator( spaces.hdiv, 1 );
    eval::SplineSpaceEvaluator l2_evaluator( spaces.l2, 1 );
    for( eval::SplineSpaceEvaluator* evaluator :
         { &geometry_evaluator,
           &h1_evaluator,
           &hcurl_evaluator,
           &hdiv_evaluator,
           &l2_evaluator } )
        evaluator->localizeElement( spaces.element );

    const Eigen::Vector2d parent_coordinates( 0.37, 0.41 );
    const Eigen::Vector2d lengths =
        spaces.h1.basisComplex()
            .parametricAtlas()
            .parametricLengths( spaces.element )
            .head<2>();

    const Eigen::MatrixXd h1_finite_difference =
        finiteDifferenceSpatialDerivatives(
            h1_evaluator,
            geometry_evaluator,
            spaces.geometry_control_points,
            parent_coordinates,
            lengths,
            [&]() { return h1_evaluator.evaluateBasisValuesAtParentPoint(); } );
    CHECK(
        ( eval::evaluateSpatialH1BasisFirstDerivatives(
              h1_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) -
          h1_finite_difference )
                .norm() < 1e-7 );

    const Eigen::MatrixXd hcurl_finite_difference =
        finiteDifferenceSpatialDerivatives(
            hcurl_evaluator,
            geometry_evaluator,
            spaces.geometry_control_points,
            parent_coordinates,
            lengths,
            [&]() {
                return eval::evaluateSpatialHCurlBasisValues(
                    hcurl_evaluator,
                    geometry_evaluator,
                    spaces.geometry_control_points );
            } );
    CHECK(
        ( eval::evaluateSpatialHCurlBasisFirstDerivatives(
              hcurl_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) -
          hcurl_finite_difference )
                .norm() < 1e-7 );

    const Eigen::MatrixXd hdiv_finite_difference =
        finiteDifferenceSpatialDerivatives(
            hdiv_evaluator,
            geometry_evaluator,
            spaces.geometry_control_points,
            parent_coordinates,
            lengths,
            [&]() {
                return eval::evaluateSpatialHDivBasisValues(
                    hdiv_evaluator,
                    geometry_evaluator,
                    spaces.geometry_control_points );
            } );
    CHECK(
        ( eval::evaluateSpatialHDivBasisFirstDerivatives(
              hdiv_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) -
          hdiv_finite_difference )
                .norm() < 1e-7 );

    const Eigen::MatrixXd l2_finite_difference =
        finiteDifferenceSpatialDerivatives(
            l2_evaluator,
            geometry_evaluator,
            spaces.geometry_control_points,
            parent_coordinates,
            lengths,
            [&]() {
                return eval::evaluateSpatialL2BasisValues(
                    l2_evaluator,
                    geometry_evaluator,
                    spaces.geometry_control_points );
            } );
    CHECK(
        ( eval::evaluateSpatialL2BasisFirstDerivatives(
              l2_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) -
          l2_finite_difference )
                .norm() < 1e-7 );
}

TEST_CASE( "Nonuniform elements preserve H(div) normal and H(curl) tangential traces" )
{
    const KnotVector knot_s( { 0, 0, 0, 1, 3, 3, 3 }, 1e-10 );
    const KnotVector knot_t( { 0, 0, 0, 2, 2, 2 }, 1e-10 );
    const TPSplineSpace h1 =
        buildBSpline( { knot_s, knot_t }, { 2, 2 } );
    const auto hdiv_complex =
        std::make_shared<const VectorConformingBasisComplex>(
            h1.basisComplexPtr(),
            ConformingType::Divergence );
    const VectorConformingTPSplineSpace hdiv( hdiv_complex, h1 );
    const auto hcurl_complex =
        std::make_shared<const VectorConformingBasisComplex>(
            h1.basisComplexPtr(),
            ConformingType::Curl );
    const VectorConformingTPSplineSpace hcurl( hcurl_complex, h1 );

    Eigen::MatrixXd geometry_control_points( 2, h1.numFunctions() );
    const std::array<double, 4> greville_s{ 0.0, 0.5, 2.0, 3.0 };
    const std::array<double, 3> greville_t{ 0.0, 1.0, 2.0 };
    for( size_t j = 0; j < greville_t.size(); ++j )
    {
        for( size_t i = 0; i < greville_s.size(); ++i )
        {
            const Eigen::Index column =
                static_cast<Eigen::Index>( j * greville_s.size() + i );
            geometry_control_points.col( column ) =
                Eigen::Vector2d( greville_s.at( i ), greville_t.at( j ) );
        }
    }

    std::vector<topology::Face> elements;
    topology::iterateCellsWhile(
        h1.basisComplex().parametricAtlas().cmap(),
        2,
        [&]( const topology::Face& face ) {
            elements.push_back( face );
            return true;
        } );
    REQUIRE( elements.size() == 2 );
    std::sort(
        elements.begin(),
        elements.end(),
        [&]( const topology::Face& first, const topology::Face& second ) {
            return h1.basisComplex()
                       .parametricAtlas()
                       .parametricStarts( first )( 0 ) <
                   h1.basisComplex()
                       .parametricAtlas()
                       .parametricStarts( second )( 0 );
        } );

    eval::SplineSpaceEvaluator geometry_evaluator( h1, 1 );
    eval::SplineSpaceEvaluator hdiv_evaluator( hdiv, 0 );
    eval::SplineSpaceEvaluator hcurl_evaluator( hcurl, 0 );
    const Eigen::Vector2d left_parent( 1.0, 0.37 );
    const Eigen::Vector2d right_parent( 0.0, 0.37 );

    const auto evaluate_side =
        [&]( const topology::Face& element,
             const Eigen::Vector2d& parent_coordinates ) {
            geometry_evaluator.localizeElement( element );
            hdiv_evaluator.localizeElement( element );
            hcurl_evaluator.localizeElement( element );
            localizeAt(
                hdiv_evaluator,
                geometry_evaluator,
                parent_coordinates );
            hcurl_evaluator.localizeParentPoint(
                parentPoint2d( parent_coordinates ) );
            return std::make_tuple(
                hdiv.connectivity( element ),
                eval::evaluateSpatialHDivBasisValues(
                    hdiv_evaluator,
                    geometry_evaluator,
                    geometry_control_points ),
                hcurl.connectivity( element ),
                eval::evaluateSpatialHCurlBasisValues(
                    hcurl_evaluator,
                    geometry_evaluator,
                    geometry_control_points ) );
        };

    const auto [left_hdiv_connectivity,
                left_hdiv_values,
                left_hcurl_connectivity,
                left_hcurl_values] =
        evaluate_side( elements.at( 0 ), left_parent );
    const auto [right_hdiv_connectivity,
                right_hdiv_values,
                right_hcurl_connectivity,
                right_hcurl_values] =
        evaluate_side( elements.at( 1 ), right_parent );

    size_t nonzero_hdiv_trace_functions = 0;
    for( size_t left = 0; left < left_hdiv_connectivity.size(); ++left )
    {
        const auto found = std::find(
            right_hdiv_connectivity.begin(),
            right_hdiv_connectivity.end(),
            left_hdiv_connectivity.at( left ) );
        if( found == right_hdiv_connectivity.end() ) continue;
        const size_t right =
            std::distance( right_hdiv_connectivity.begin(), found );
        CHECK(
            std::abs(
                left_hdiv_values( left, 0 ) -
                right_hdiv_values( right, 0 ) ) < 1e-12 );
        if( std::abs( left_hdiv_values( left, 0 ) ) > 1e-12 )
            ++nonzero_hdiv_trace_functions;
    }
    CHECK( nonzero_hdiv_trace_functions > 0 );

    size_t nonzero_hcurl_trace_functions = 0;
    for( size_t left = 0; left < left_hcurl_connectivity.size(); ++left )
    {
        const auto found = std::find(
            right_hcurl_connectivity.begin(),
            right_hcurl_connectivity.end(),
            left_hcurl_connectivity.at( left ) );
        if( found == right_hcurl_connectivity.end() ) continue;
        const size_t right =
            std::distance( right_hcurl_connectivity.begin(), found );
        CHECK(
            std::abs(
                left_hcurl_values( left, 1 ) -
                right_hcurl_values( right, 1 ) ) < 1e-12 );
        if( std::abs( left_hcurl_values( left, 1 ) ) > 1e-12 )
            ++nonzero_hcurl_trace_functions;
    }
    CHECK( nonzero_hcurl_trace_functions > 0 );
}

namespace
{
    param::ParentPoint parentPoint3d( const Eigen::Vector3d& coordinates )
    {
        return param::ParentPoint(
            param::cubeDomain( 3 ),
            coordinates,
            { false, false, false, false, false, false } );
    }

    void localizeAt3d(
        eval::SplineSpaceEvaluator& basis_evaluator,
        eval::SplineSpaceEvaluator& geometry_evaluator,
        const Eigen::Vector3d& coordinates )
    {
        const param::ParentPoint point = parentPoint3d( coordinates );
        basis_evaluator.localizeParentPoint( point );
        geometry_evaluator.localizeParentPoint( point );
    }

    template <typename ValueEvaluator>
    Eigen::MatrixXd finiteDifferenceSpatialDerivatives3d(
        eval::SplineSpaceEvaluator& basis_evaluator,
        eval::SplineSpaceEvaluator& geometry_evaluator,
        const Eigen::MatrixXd& geometry_control_points,
        const Eigen::Vector3d& parent_coordinates,
        const Eigen::Vector3d& parametric_lengths,
        ValueEvaluator&& evaluate_values )
    {
        constexpr double epsilon = 1e-6;
        std::array<Eigen::MatrixXd, 3> parametric_derivatives;
        for( size_t direction = 0; direction < 3; ++direction )
        {
            Eigen::Vector3d parent_step = Eigen::Vector3d::Zero();
            parent_step( direction ) =
                epsilon / parametric_lengths( direction );
            localizeAt3d(
                basis_evaluator,
                geometry_evaluator,
                parent_coordinates + parent_step );
            const Eigen::MatrixXd plus = evaluate_values();
            localizeAt3d(
                basis_evaluator,
                geometry_evaluator,
                parent_coordinates - parent_step );
            const Eigen::MatrixXd minus = evaluate_values();
            parametric_derivatives.at( direction ) =
                ( plus - minus ) / ( 2.0 * epsilon );
        }

        localizeAt3d(
            basis_evaluator,
            geometry_evaluator,
            parent_coordinates );
        const Eigen::MatrixXd inverse_jacobian =
            geometry_evaluator
                .evaluateParametricToSpatialJacobian(
                    geometry_control_points )
                .inverse();
        const Eigen::Index n_functions =
            parametric_derivatives.front().rows();
        const Eigen::Index n_components =
            parametric_derivatives.front().cols();
        Eigen::MatrixXd output(
            n_functions,
            n_components * inverse_jacobian.cols() );
        for( Eigen::Index function = 0; function < n_functions; ++function )
        {
            Eigen::MatrixXd derivative_wrt_parametric( n_components, 3 );
            for( size_t direction = 0; direction < 3; ++direction )
                derivative_wrt_parametric.col( direction ) =
                    parametric_derivatives.at( direction )
                        .row( function )
                        .transpose();
            const Eigen::MatrixXd derivative_wrt_spatial =
                derivative_wrt_parametric * inverse_jacobian;
            for( Eigen::Index direction = 0;
                 direction < derivative_wrt_spatial.cols();
                 ++direction )
            {
                for( Eigen::Index component = 0;
                     component < derivative_wrt_spatial.rows();
                     ++component )
                {
                    output(
                        function,
                        direction * n_components + component ) =
                        derivative_wrt_spatial( component, direction );
                }
            }
        }
        return output;
    }

    struct SpatialTestSpaces3d
    {
        SpatialTestSpaces3d()
            : h1( buildBSpline(
                  { KnotVector( { 0, 0, 0, 2, 2, 2 }, 1e-10 ),
                    KnotVector( { 0, 0, 0, 3, 3, 3 }, 1e-10 ),
                    KnotVector( { 0, 0, 0, 4, 4, 4 }, 1e-10 ) },
                  { 2, 2, 2 } ) ),
              hcurl_complex( std::make_shared<const VectorConformingBasisComplex>(
                  h1.basisComplexPtr(),
                  ConformingType::Curl ) ),
              hcurl( hcurl_complex, h1 ),
              hdiv_complex( std::make_shared<const VectorConformingBasisComplex>(
                  h1.basisComplexPtr(),
                  ConformingType::Divergence ) ),
              hdiv( hdiv_complex, h1 ),
              geometry_control_points( 3, h1.numFunctions() )
        {
            for( size_t k = 0; k < 3; ++k )
            {
                for( size_t j = 0; j < 3; ++j )
                {
                    for( size_t i = 0; i < 3; ++i )
                    {
                        const double s_control =
                            0.5 * static_cast<double>( i );
                        const double t_control =
                            0.5 * static_cast<double>( j );
                        const double u_control =
                            0.5 * static_cast<double>( k );
                        const Eigen::Index column =
                            static_cast<Eigen::Index>(
                                9 * k + 3 * j + i );
                        geometry_control_points( 0, column ) =
                            s_control +
                            0.1 * s_control * t_control;
                        geometry_control_points( 1, column ) =
                            t_control +
                            0.08 * t_control * u_control;
                        geometry_control_points( 2, column ) =
                            u_control + ( i == 2 ? 0.06 : 0.0 );
                    }
                }
            }
        }

        TPSplineSpace h1;
        std::shared_ptr<const VectorConformingBasisComplex> hcurl_complex;
        VectorConformingTPSplineSpace hcurl;
        std::shared_ptr<const VectorConformingBasisComplex> hdiv_complex;
        VectorConformingTPSplineSpace hdiv;
        Eigen::MatrixXd geometry_control_points;
        topology::Volume element{ topology::Dart( 0 ) };
    };
}

TEST_CASE( "Three-dimensional spatial derivatives agree with finite differences" )
{
    SpatialTestSpaces3d spaces;
    eval::SplineSpaceEvaluator geometry_evaluator( spaces.h1, 2 );
    eval::SplineSpaceEvaluator h1_evaluator( spaces.h1, 1 );
    eval::SplineSpaceEvaluator hcurl_evaluator( spaces.hcurl, 1 );
    eval::SplineSpaceEvaluator hdiv_evaluator( spaces.hdiv, 1 );
    for( eval::SplineSpaceEvaluator* evaluator :
         { &geometry_evaluator,
           &h1_evaluator,
           &hcurl_evaluator,
           &hdiv_evaluator } )
        evaluator->localizeElement( spaces.element );

    const Eigen::Vector3d parent_coordinates( 0.37, 0.41, 0.53 );
    const Eigen::Vector3d lengths =
        spaces.h1.basisComplex()
            .parametricAtlas()
            .parametricLengths( spaces.element )
            .head<3>();
    localizeAt3d(
        h1_evaluator,
        geometry_evaluator,
        parent_coordinates );

    const std::vector<Eigen::MatrixXd> analytic_jacobian_derivatives =
        eval::evaluateParametricToSpatialJacobianFirstDerivatives(
            geometry_evaluator,
            spaces.geometry_control_points );
    const Eigen::VectorXd analytic_determinant_derivatives =
        eval::evaluateParametricToSpatialJacobianDeterminantFirstDerivatives(
            geometry_evaluator,
            spaces.geometry_control_points );
    constexpr double epsilon = 1e-6;
    for( size_t direction = 0; direction < 3; ++direction )
    {
        Eigen::Vector3d parent_step = Eigen::Vector3d::Zero();
        parent_step( direction ) = epsilon / lengths( direction );
        geometry_evaluator.localizeParentPoint(
            parentPoint3d( parent_coordinates + parent_step ) );
        const Eigen::MatrixXd plus =
            geometry_evaluator.evaluateParametricToSpatialJacobian(
                spaces.geometry_control_points );
        geometry_evaluator.localizeParentPoint(
            parentPoint3d( parent_coordinates - parent_step ) );
        const Eigen::MatrixXd minus =
            geometry_evaluator.evaluateParametricToSpatialJacobian(
                spaces.geometry_control_points );
        CHECK(
            ( analytic_jacobian_derivatives.at( direction ) -
              ( plus - minus ) / ( 2.0 * epsilon ) )
                    .norm() < 1e-8 );
        CHECK(
            std::abs(
                analytic_determinant_derivatives( direction ) -
                ( plus.determinant() - minus.determinant() ) /
                    ( 2.0 * epsilon ) ) < 1e-8 );
    }

    const Eigen::MatrixXd h1_finite_difference =
        finiteDifferenceSpatialDerivatives3d(
            h1_evaluator,
            geometry_evaluator,
            spaces.geometry_control_points,
            parent_coordinates,
            lengths,
            [&]() { return h1_evaluator.evaluateBasisValuesAtParentPoint(); } );
    CHECK(
        ( eval::evaluateSpatialH1BasisFirstDerivatives(
              h1_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) -
          h1_finite_difference )
                .norm() < 2e-7 );

    const Eigen::MatrixXd hcurl_finite_difference =
        finiteDifferenceSpatialDerivatives3d(
            hcurl_evaluator,
            geometry_evaluator,
            spaces.geometry_control_points,
            parent_coordinates,
            lengths,
            [&]() {
                return eval::evaluateSpatialHCurlBasisValues(
                    hcurl_evaluator,
                    geometry_evaluator,
                    spaces.geometry_control_points );
            } );
    CHECK(
        ( eval::evaluateSpatialHCurlBasisFirstDerivatives(
              hcurl_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) -
          hcurl_finite_difference )
                .norm() < 2e-7 );

    const Eigen::MatrixXd hdiv_finite_difference =
        finiteDifferenceSpatialDerivatives3d(
            hdiv_evaluator,
            geometry_evaluator,
            spaces.geometry_control_points,
            parent_coordinates,
            lengths,
            [&]() {
                return eval::evaluateSpatialHDivBasisValues(
                    hdiv_evaluator,
                    geometry_evaluator,
                    spaces.geometry_control_points );
            } );
    CHECK(
        ( eval::evaluateSpatialHDivBasisFirstDerivatives(
              hdiv_evaluator,
              geometry_evaluator,
              spaces.geometry_control_points ) -
          hdiv_finite_difference )
                .norm() < 2e-7 );
}

Eigen::MatrixXd quadballControlPoints()
{
    // See table 4 of "Tiling the Sphere with Rational Bezier Patches" (James E. Cobb, 1988)
    // Available here: https://collections.lib.utah.edu/dl_files/4e/77/4e7746dd53c79f8557272b92b47d2d407da4931a.pdf

    const double rt2 = std::sqrt(2.0);
    const double rt3 = std::sqrt(3.0);
    const double rt6 = std::sqrt(6.0);

    Eigen::MatrixXd P(4,25);

    P.row(0) = ( Eigen::RowVectorXd(25) <<
        4*(1-rt3),
        -rt2,
        0,
        rt2,
        4*(rt3-1),

        rt2*(rt3-4),
        (2-3*rt3)/2.0,
        0,
        (3*rt3-2)/2.0,
        rt2*(4-rt3),

        4*(1-2*rt3)/3.0,
        rt2*(2*rt3-7)/3.0,
        0,
        rt2*(7-2*rt3)/3.0,
        4*(2*rt3-1)/3.0,

        rt2*(rt3-4),
        (2-3*rt3)/2.0,
        0,
        (3*rt3-2)/2.0,
        rt2*(4-rt3),

        4*(1-rt3),
        -rt2,
        0,
        rt2,
        4*(rt3-1)
    ).finished();

    P.row(1) = ( Eigen::RowVectorXd(25) <<
        4*(1-rt3),
        rt2*(rt3-4),
        4*(1-2*rt3)/3.0,
        rt2*(rt3-4),
        4*(1-rt3),

        -rt2,
        (2-3*rt3)/2.0,
        rt2*(2*rt3-7)/3.0,
        (2-3*rt3)/2.0,
        -rt2,

        0,
        0,
        0,
        0,
        0,

        rt2,
        (3*rt3-2)/2.0,
        rt2*(7-2*rt3)/3.0,
        (3*rt3-2)/2.0,
        rt2,

        4*(rt3-1),
        rt2*(4-rt3),
        4*(2*rt3-1)/3.0,
        rt2*(4-rt3),
        4*(rt3-1)
    ).finished();

    P.row(2) = ( Eigen::RowVectorXd(25) <<
        4*(1-rt3),
        rt2*(rt3-4),
        4*(1-2*rt3)/3.0,
        rt2*(rt3-4),
        4*(1-rt3),

        rt2*(rt3-4),
        -(rt3+6)/2.0,
        -5*rt6/3.0,
        -(rt3+6)/2.0,
        rt2*(rt3-4),

        4*(1-2*rt3)/3.0,
        -5*rt6/3.0,
        4*(rt3-5)/3.0,
        -5*rt6/3.0,
        4*(1-2*rt3)/3.0,

        rt2*(rt3-4),
        -(rt3+6)/2.0,
        -5*rt6/3.0,
        -(rt3+6)/2.0,
        rt2*(rt3-4),

        4*(1-rt3),
        rt2*(rt3-4),
        4*(1-2*rt3)/3.0,
        rt2*(rt3-4),
        4*(1-rt3)
    ).finished();

    P.row(3) = ( Eigen::RowVectorXd(25) <<
        4*(3-rt3),
        rt2*(3*rt3-2),
        4*(5-rt3)/3.0,
        rt2*(3*rt3-2),
        4*(3-rt3),

        rt2*(3*rt3-2),
        (rt3+6)/2.0,
        rt2*(rt3+6)/3.0,
        (rt3+6)/2.0,
        rt2*(3*rt3-2),

        4*(5-rt3)/3.0,
        rt2*(rt3+6)/3.0,
        4*(5*rt3-1)/9.0,
        rt2*(rt3+6)/3.0,
        4*(5-rt3)/3.0,

        rt2*(3*rt3-2),
        (rt3+6)/2.0,
        rt2*(rt3+6)/3.0,
        (rt3+6)/2.0,
        rt2*(3*rt3-2),

        4*(3-rt3),
        rt2*(3*rt3-2),
        4*(5-rt3)/3.0,
        rt2*(3*rt3-2),
        4*(3-rt3)
    ).finished();

    return P;
}


TEST_CASE( "NURBS Spline space evaluation" )
{
    const KnotVector kv( {0,0,0,0,0,1,1,1,1,1}, 1e-10 );
    const size_t degree = 4;

    const auto cmap = std::make_shared<const CombinatorialMap1d>( numElements( kv ) );
    const auto param = std::make_shared<const ParametricAtlas1d>( cmap, parametricLengths( kv ) );
    const auto bc = std::make_shared<const BasisComplex1d>( param, degree );
    const auto ss = std::make_shared<const BSplineSpace1d>( bc, kv );

    const auto cmap_2d = std::make_shared<const TPCombinatorialMap>( cmap, cmap );
    const auto param_2d = std::make_shared<const TPParametricAtlas>( cmap_2d, param, param );
    const auto bc_2d = std::make_shared<const TPBasisComplex>( param_2d, bc, bc );
    const TPSplineSpace ss_2d( bc_2d, ss, ss );

    const Eigen::MatrixXd control_points = quadballControlPoints();

    eval::NURBSSpaceEvaluator nurbs_evals( ss_2d , 2 );

    iterateCellsWhile( *cmap_2d, 2, [&]( const topology::Face& f ) {
        const ParentPoint ppt( param_2d->parentDomain( f ), Eigen::Vector2d(0.5,0.5), BaryCoordIsZeroVec{false,false,false,false} );
        nurbs_evals.localizeElement( f );
        nurbs_evals.localizeParentPoint( ppt );
        const Eigen::Vector3d eval = nurbs_evals.evaluateManifold( control_points );
        const Eigen::MatrixXd jac = nurbs_evals.evaluateParentToSpatialJacobian( control_points );
        const Eigen::MatrixXd hess = nurbs_evals.evaluateParentToSpatialHessian( control_points );

        LOG( LOG_EVALS ) << "Evaluation at face " << f.dart() << ":\n";
        LOG( LOG_EVALS ) << "  Point: " << eval.transpose() << "\n";
        LOG( LOG_EVALS ) << "  Jacobian:\n" << jac << "\n";
        LOG( LOG_EVALS ) << "  Hessian:\n" << hess << "\n";

        CHECK( util::equals( eval.norm(), 1.0, 1e-14 ) );
        CHECK( util::equals( jac.col(0).norm(), jac.col(1).norm(), 1e-14 ) );
        CHECK( util::equals( jac.col(0).norm(), jac(0, 0), 1e-14 ) );
        CHECK( util::equals( jac.col(1).norm(), jac(1, 1), 1e-14 ) );
        CHECK( util::equals( jac.col(0).dot( jac.col(1) ), 0.0, 1e-14 ) );

        // I don't have good intuition for why this should be true nor the time to check it,
        // but it is true, so we'll check it just for the sake of regression prevention.
        CHECK( util::equals(hess.col(0).norm(), hess.col(2).norm(), 1e-14 ) );
        CHECK( util::equals(hess.col(1).norm(), 0, 1e-14 ) );
        
        return true;
    } );
}
