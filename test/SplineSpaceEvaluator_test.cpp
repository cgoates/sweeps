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
            primal_evals.localizePoint( ppt );

            const auto conn = ss_2d.connectivity( f );
            const auto basis = primal_evals.evaluateBasis();
            const auto derv = primal_evals.evaluateFirstDerivatives();

            const ParentPoint ppt_opp =
                pointOnBoundary( param_2d->parentDomain( f_opp ),
                                 parentDomainBoundary( *param_2d, topology::Edge( maybe_opp_e.value() ) ) );
            primal_evals.localizeElement( f_opp );
            primal_evals.localizePoint( ppt_opp );

            const auto conn_opp = ss_2d.connectivity( f_opp );
            const auto basis_opp = primal_evals.evaluateBasis();
            const auto derv_opp = primal_evals.evaluateFirstDerivatives();

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

TEST_CASE( "Piola transforms are invariant to parent-to-parametric interval scaling" )
{
    const double ptol = 1e-10;
    const size_t degree = 1;

    const auto build_discretization = [&]( const KnotVector& kv_s, const KnotVector& kv_t ) {
        const TPSplineSpace h1 = buildBSpline( { kv_s, kv_t }, { degree, degree } );
        const auto hdiv_bc = std::make_shared<const VectorConformingBasisComplex>( h1.basisComplexPtr() );
        const VectorConformingTPSplineSpace hdiv( hdiv_bc, h1 );
        const auto l2_bc =
            std::make_shared<const TPBasisComplex>( h1.basisComplex().parametricAtlasPtr(),
                                                   hdiv.reducedDegree1dBases().at( 0 )->basisComplexPtr(),
                                                   hdiv.reducedDegree1dBases().at( 1 )->basisComplexPtr() );
        const TPSplineSpace l2( l2_bc, hdiv.reducedDegree1dBases().at( 0 ), hdiv.reducedDegree1dBases().at( 1 ) );
        return std::tuple<TPSplineSpace, VectorConformingTPSplineSpace, TPSplineSpace>( h1, hdiv, l2 );
    };

    const KnotVector kv_unit( { 0, 0, 1, 1 }, ptol );
    const KnotVector kv_double( { 0, 0, 2, 2 }, ptol );

    auto [h1_unit, hdiv_unit, l2_unit] = build_discretization( kv_unit, kv_unit );
    auto [h1_scaled, hdiv_scaled, l2_scaled] = build_discretization( kv_double, kv_unit );

    eval::SplineSpaceEvaluator h1_eval_unit( h1_unit, 2 );
    eval::SplineSpaceEvaluator hdiv_eval_unit( hdiv_unit, 1 );
    eval::SplineSpaceEvaluator l2_eval_unit( l2_unit, 1 );

    eval::SplineSpaceEvaluator h1_eval_scaled( h1_scaled, 2 );
    eval::SplineSpaceEvaluator hdiv_eval_scaled( hdiv_scaled, 1 );
    eval::SplineSpaceEvaluator l2_eval_scaled( l2_scaled, 1 );

    const topology::Face elem( topology::Dart( 0 ) );
    const param::ParentPoint pt( param::cubeDomain( 2 ), Eigen::Vector2d( 0.37, 0.61 ), { false, false, false, false } );

    h1_eval_unit.localizeElement( elem );
    hdiv_eval_unit.localizeElement( elem );
    l2_eval_unit.localizeElement( elem );
    h1_eval_scaled.localizeElement( elem );
    hdiv_eval_scaled.localizeElement( elem );
    l2_eval_scaled.localizeElement( elem );

    h1_eval_unit.localizePoint( pt );
    hdiv_eval_unit.localizePoint( pt );
    l2_eval_unit.localizePoint( pt );
    h1_eval_scaled.localizePoint( pt );
    hdiv_eval_scaled.localizePoint( pt );
    l2_eval_scaled.localizePoint( pt );

    Eigen::MatrixX2d unit_geom( 4, 2 );
    unit_geom << 0.0, 0.0,
                 1.0, 0.0,
                 0.0, 1.0,
                 1.0, 1.0;
    const Eigen::MatrixX2d scaled_geom = unit_geom;

    CHECK( ( eval::piolaTransformedH1FirstDerivatives( h1_eval_unit, h1_eval_unit, unit_geom.transpose() ) -
             eval::piolaTransformedH1FirstDerivatives( h1_eval_scaled, h1_eval_scaled, scaled_geom.transpose() ) )
               .norm() < 1e-12 );
    CHECK( ( eval::piolaTransformedHDivBasis( hdiv_eval_unit, h1_eval_unit, unit_geom.transpose() ) -
             eval::piolaTransformedHDivBasis( hdiv_eval_scaled, h1_eval_scaled, scaled_geom.transpose() ) )
               .norm() < 1e-12 );
    CHECK( ( eval::piolaTransformedHDivFirstDerivatives( hdiv_eval_unit, h1_eval_unit, unit_geom.transpose() ) -
             eval::piolaTransformedHDivFirstDerivatives( hdiv_eval_scaled, h1_eval_scaled, scaled_geom.transpose() ) )
               .norm() < 1e-12 );
    CHECK( ( eval::piolaTransformedL2Basis( l2_eval_unit, h1_eval_unit, unit_geom.transpose() ) -
             eval::piolaTransformedL2Basis( l2_eval_scaled, h1_eval_scaled, scaled_geom.transpose() ) )
               .norm() < 1e-12 );
    CHECK( ( eval::piolaTransformedL2FirstDerivatives( l2_eval_unit, h1_eval_unit, unit_geom.transpose() ) -
             eval::piolaTransformedL2FirstDerivatives( l2_eval_scaled, h1_eval_scaled, scaled_geom.transpose() ) )
               .norm() < 1e-12 );
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
        nurbs_evals.localizePoint( ppt );
        const Eigen::Vector3d eval = nurbs_evals.evaluateManifold( control_points );
        const Eigen::MatrixXd jac = nurbs_evals.evaluateJacobian( control_points );
        const Eigen::MatrixXd hess = nurbs_evals.evaluateHessian( control_points );

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
