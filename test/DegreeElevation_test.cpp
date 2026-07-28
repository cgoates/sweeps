#include <catch2/catch_test_macros.hpp>
#include <BezierTransfer.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParentDomain.hpp>
#include <ParentPoint.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <TPSplineSpace.hpp>
#include <cmath>

using namespace basis;
using namespace param;
using namespace topology;

namespace
{
    long double logBinomial( const size_t n, const size_t k )
    {
        return std::lgamma( static_cast<long double>( n ) + 1.0L ) -
               std::lgamma( static_cast<long double>( k ) + 1.0L ) -
               std::lgamma( static_cast<long double>( n - k ) + 1.0L );
    }

    Eigen::VectorXd bernsteinValues( const size_t degree, const double x )
    {
        Eigen::VectorXd out( degree + 1 );
        for( size_t i = 0; i <= degree; i++ )
        {
            const long double log_value =
                logBinomial( degree, i ) +
                static_cast<long double>( i ) * std::log( static_cast<long double>( x ) ) +
                static_cast<long double>( degree - i ) * std::log( static_cast<long double>( 1.0 - x ) );
            out( i ) = static_cast<double>( std::exp( log_value ) );
        }
        return out;
    }

    Eigen::MatrixXd deterministicControlPoints( const size_t spatial_dim, const size_t num_functions )
    {
        Eigen::MatrixXd cpts( spatial_dim, num_functions );
        for( size_t i = 0; i < spatial_dim; i++ )
        {
            for( size_t j = 0; j < num_functions; j++ )
            {
                cpts( i, j ) = std::sin( 0.37 * static_cast<double>( i + 1 ) * static_cast<double>( j + 1 ) ) +
                               0.11 * static_cast<double>( i ) +
                               0.03 * static_cast<double>( j );
            }
        }
        return cpts;
    }

    void checkLocalDegreeElevationIdentity( const KnotVector& kv,
                                            const size_t p,
                                            const size_t q,
                                            const double ptol,
                                            const double check_tol )
    {
        const Eigen::SparseMatrix<double> T = degreeElevationOp( kv, p, q, ptol );
        const KnotVector elevated_kv = degreeElevatedKnotVector( kv, q - p );
        const std::vector<LocalExtraction> source_locals = localExtractions( kv, p );
        const std::vector<LocalExtraction> target_locals = localExtractions( elevated_kv, q );
        const Eigen::MatrixXd M = bernsteinDegreeElevationMatrix( p, q );

        REQUIRE( source_locals.size() == target_locals.size() );
        for( size_t elem_ii = 0; elem_ii < source_locals.size(); elem_ii++ )
        {
            const LocalExtraction& source_local = source_locals.at( elem_ii );
            const LocalExtraction& target_local = target_locals.at( elem_ii );
            Eigen::MatrixXd local_T( source_local.connectivity.size(), target_local.connectivity.size() );

            for( Eigen::Index i = 0; i < local_T.rows(); i++ )
            {
                for( Eigen::Index j = 0; j < local_T.cols(); j++ )
                {
                    local_T( i, j ) =
                        T.coeff( source_local.connectivity.at( i ), target_local.connectivity.at( j ) );
                }
            }

            CHECK( ( local_T * target_local.extraction - source_local.extraction * M ).norm() < check_tol );
        }
    }
}

TEST_CASE( "Bernstein degree elevation matrix preserves Bernstein values" )
{
    const size_t p = 30;
    const size_t q = 45;
    const Eigen::MatrixXd M = bernsteinDegreeElevationMatrix( p, q );

    CHECK( M.rows() == p + 1 );
    CHECK( M.cols() == q + 1 );

    for( Eigen::Index i = 0; i < M.rows(); i++ )
    {
        for( Eigen::Index j = 0; j < M.cols(); j++ )
        {
            CHECK( std::isfinite( M( i, j ) ) );
            CHECK( M( i, j ) >= -1e-15 );
        }
    }

    for( const double x : { 0.07, 0.23, 0.51, 0.89 } )
    {
        const Eigen::VectorXd elevated_values = M * bernsteinValues( q, x );
        CHECK( ( elevated_values - bernsteinValues( p, x ) ).norm() < 1e-12 );
    }
}

TEST_CASE( "Univariate degree elevation operator satisfies local extraction identity" )
{
    const double ptol = 1e-10;
    const KnotVector kv( { 0, 0, 0, 0, 0.5, 1, 2, 2, 2, 2 }, ptol );

    checkLocalDegreeElevationIdentity( kv, 3, 5, ptol, 1e-10 );
}

TEST_CASE( "High-degree univariate degree elevation operator satisfies local extraction identity" )
{
    const double ptol = 1e-10;
    const KnotVector kv( { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0.25, 0.5, 0.5, 0.75,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                         ptol );

    checkLocalDegreeElevationIdentity( kv, 12, 20, ptol, 1e-8 );
}

TEST_CASE( "Tensor-product degree elevation preserves surface geometry and parent derivatives" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s( { 0, 0, 0, 1, 2, 2, 2 }, ptol );
    const KnotVector kv_t( { 0, 0, 0, 0, 0.5, 1, 1, 1, 1 }, ptol );
    const SmallVector<size_t, 3> source_degrees( { 2, 3 } );
    const SmallVector<size_t, 3> target_degrees( { 4, 4 } );

    const TPSplineSpace source = buildBSpline( { kv_s, kv_t }, source_degrees );
    const TPSplineSpace target = degreeElevatedSpace( source, target_degrees );
    const Eigen::SparseMatrix<double> T = degreeElevationOp( source, target_degrees, ptol );

    const Eigen::MatrixXd source_cpts = deterministicControlPoints( 3, source.numFunctions() );
    const Eigen::MatrixXd target_cpts = source_cpts * T;

    eval::SplineSpaceEvaluator source_evaler( source, 1 );
    eval::SplineSpaceEvaluator target_evaler( target, 1 );
    const ParentPoint ppt( cubeDomain( 2 ), Eigen::Vector2d( 0.37, 0.61 ), { false, false, false, false } );

    iterateCellsWhile( source.basisComplex().parametricAtlas().cmap(), 2, [&]( const Cell& c ) {
        source_evaler.localizeElement( c );
        target_evaler.localizeElement( c );
        source_evaler.localizeParentPoint( ppt );
        target_evaler.localizeParentPoint( ppt );

        CHECK( ( source_evaler.evaluateManifold( source_cpts ) -
                 target_evaler.evaluateManifold( target_cpts ) ).norm() < 1e-10 );
        CHECK( ( source_evaler.evaluateParentToSpatialJacobian( source_cpts ) -
                 target_evaler.evaluateParentToSpatialJacobian( target_cpts ) ).norm() < 1e-10 );
        return true;
    } );
}

TEST_CASE( "Tensor-product degree elevation preserves volume geometry" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s( { 0, 0, 1, 2, 2 }, ptol );
    const KnotVector kv_t( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector kv_u( { 0, 0, 1, 1 }, ptol );
    const SmallVector<size_t, 3> source_degrees( { 1, 2, 1 } );
    const SmallVector<size_t, 3> target_degrees( { 2, 3, 3 } );

    const TPSplineSpace source = buildBSpline( { kv_s, kv_t, kv_u }, source_degrees );
    const TPSplineSpace target = degreeElevatedSpace( source, target_degrees );
    const Eigen::SparseMatrix<double> T = degreeElevationOp( source, target_degrees, ptol );

    const Eigen::MatrixXd source_cpts = deterministicControlPoints( 3, source.numFunctions() );
    const Eigen::MatrixXd target_cpts = source_cpts * T;

    eval::SplineSpaceEvaluator source_evaler( source, 0 );
    eval::SplineSpaceEvaluator target_evaler( target, 0 );
    const ParentPoint ppt(
        cubeDomain( 3 ), Eigen::Vector3d( 0.25, 0.5, 0.75 ), { false, false, false, false, false, false } );

    iterateCellsWhile( source.basisComplex().parametricAtlas().cmap(), 3, [&]( const Cell& c ) {
        source_evaler.localizeElement( c );
        target_evaler.localizeElement( c );
        source_evaler.localizeParentPoint( ppt );
        target_evaler.localizeParentPoint( ppt );

        CHECK( ( source_evaler.evaluateManifold( source_cpts ) -
                 target_evaler.evaluateManifold( target_cpts ) ).norm() < 1e-10 );
        return true;
    } );
}
