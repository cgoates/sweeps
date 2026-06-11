#include <catch2/catch_test_macros.hpp>
#include <BezierTransfer.hpp>
#include <BSplineSpace1d.hpp>
#include <BasisComplex1d.hpp>
#include <CombinatorialMap1d.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParametricAtlas1d.hpp>
#include <TPSplineSpace.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

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

    double bsplineBasisValue( const KnotVector& kv, const size_t degree, const size_t function_ii, const double x )
    {
        if( degree == 0 )
            return kv.knot( function_ii ) <= x and x < kv.knot( function_ii + 1 ) ? 1.0 : 0.0;

        const double left_denominator = kv.knot( function_ii + degree ) - kv.knot( function_ii );
        const double right_denominator = kv.knot( function_ii + degree + 1 ) - kv.knot( function_ii + 1 );

        const double left =
            left_denominator == 0.0
                ? 0.0
                : ( x - kv.knot( function_ii ) ) / left_denominator *
                      bsplineBasisValue( kv, degree - 1, function_ii, x );
        const double right =
            right_denominator == 0.0
                ? 0.0
                : ( kv.knot( function_ii + degree + 1 ) - x ) / right_denominator *
                      bsplineBasisValue( kv, degree - 1, function_ii + 1, x );

        return left + right;
    }

    std::set<size_t> idSet( const std::vector<FunctionId>& conn )
    {
        std::set<size_t> out;
        std::transform( conn.begin(), conn.end(), std::inserter( out, out.end() ), []( const FunctionId& fid ) {
            return fid.id();
        } );
        return out;
    }

    size_t intersectionSize( const std::vector<FunctionId>& a, const std::vector<FunctionId>& b )
    {
        const std::set<size_t> a_ids = idSet( a );
        const std::set<size_t> b_ids = idSet( b );
        std::vector<size_t> intersection;
        std::set_intersection(
            a_ids.begin(), a_ids.end(), b_ids.begin(), b_ids.end(), std::back_inserter( intersection ) );
        return intersection.size();
    }

    void checkExtractionValues( const KnotVector& kv,
                                const size_t degree,
                                const size_t elem_ii,
                                const std::vector<FunctionId>& connectivity,
                                const Eigen::MatrixXd& extraction,
                                const double tol )
    {
        const std::vector<double> unique_knots = kv.uniqueKnots();

        REQUIRE( connectivity.size() == degree + 1 );
        REQUIRE( extraction.rows() == static_cast<Eigen::Index>( degree + 1 ) );
        REQUIRE( extraction.cols() == static_cast<Eigen::Index>( degree + 1 ) );

        for( const double parent_x : { 0.17, 0.49, 0.83 } )
        {
            const double parametric_x =
                unique_knots.at( elem_ii ) +
                parent_x * ( unique_knots.at( elem_ii + 1 ) - unique_knots.at( elem_ii ) );
            Eigen::VectorXd expected( connectivity.size() );
            for( Eigen::Index i = 0; i < expected.size(); i++ )
            {
                expected( i ) =
                    bsplineBasisValue( kv, degree, connectivity.at( static_cast<size_t>( i ) ).id(), parametric_x );
            }

            CHECK( ( extraction * bernsteinValues( degree, parent_x ) - expected ).norm() < tol );
        }
    }

    void checkLocalExtractionSharing( const KnotVector& kv,
                                      const size_t degree,
                                      const size_t expected_adjacent_overlap,
                                      const bool expect_identity_extraction )
    {
        const std::vector<LocalExtraction> locals = localExtractions( kv, degree );
        REQUIRE( locals.size() == numElements( kv ) );

        for( const LocalExtraction& local : locals )
        {
            REQUIRE( local.connectivity.size() == degree + 1 );
            CHECK( ( local.extraction.colwise().sum().transpose() -
                     Eigen::VectorXd::Ones( local.extraction.cols() ) ).norm() < 1e-14 );

            if( expect_identity_extraction )
                CHECK( ( local.extraction - Eigen::MatrixXd::Identity( degree + 1, degree + 1 ) ).norm() < 1e-14 );
        }

        for( size_t elem_ii = 0; elem_ii + 1 < locals.size(); elem_ii++ )
        {
            CHECK( intersectionSize( locals.at( elem_ii ).connectivity, locals.at( elem_ii + 1 ).connectivity ) ==
                   expected_adjacent_overlap );
        }
    }
}

TEST_CASE( "Element Bezier column offsets follow global extraction layout" )
{
    const double ptol = 1e-10;

    CHECK( elementBezierColumnOffsets( KnotVector( { 0, 0, 0, 1, 2, 3, 3, 3 }, ptol ), 2 ) ==
           std::vector<size_t>{ 0, 2, 4 } );
    CHECK( elementBezierColumnOffsets( KnotVector( { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 }, ptol ), 2 ) ==
           std::vector<size_t>{ 0, 2, 4 } );
    CHECK( elementBezierColumnOffsets( KnotVector( { 0, 0, 1, 1, 2, 2, 3, 3 }, ptol ), 1 ) ==
           std::vector<size_t>{ 0, 2, 4 } );
    CHECK( elementBezierColumnOffsets( KnotVector( { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3 }, ptol ), 2 ) ==
           std::vector<size_t>{ 0, 3, 6 } );
    CHECK( elementBezierColumnOffsets( KnotVector( { 0, 0, 0, 1, 2, 2, 3, 3, 3, 4, 4, 4 }, ptol ), 2 ) ==
           std::vector<size_t>{ 0, 2, 4, 7 } );
}

TEST_CASE( "Local extractions respect smooth, C0, and broken 1D continuity" )
{
    const double ptol = 1e-10;

    SECTION( "smooth quadratic splines share degree DOFs between adjacent elements" )
    {
        const KnotVector kv( { 0, 0, 0, 1, 2, 3, 3, 3 }, ptol );
        checkLocalExtractionSharing( kv, 2, 2, false );
    }

    SECTION( "C0 quadratic splines share exactly one DOF between adjacent elements" )
    {
        const KnotVector kv( { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 }, ptol );
        checkLocalExtractionSharing( kv, 2, 1, true );
    }

    SECTION( "broken linear splines share no DOFs between adjacent elements" )
    {
        const KnotVector kv( { 0, 0, 1, 1, 2, 2, 3, 3 }, ptol );
        checkLocalExtractionSharing( kv, 1, 0, true );
    }

    SECTION( "broken quadratic splines share no DOFs between adjacent elements" )
    {
        const KnotVector kv( { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3 }, ptol );
        checkLocalExtractionSharing( kv, 2, 0, true );
    }

    SECTION( "mixed quadratic knot multiplicities use side-specific Bezier columns" )
    {
        const size_t degree = 2;
        const KnotVector kv( { 0, 0, 0, 1, 2, 2, 3, 3, 3, 4, 4, 4 }, ptol );
        const std::vector<LocalExtraction> locals = localExtractions( kv, degree );

        REQUIRE( locals.size() == 4 );
        CHECK( intersectionSize( locals.at( 0 ).connectivity, locals.at( 1 ).connectivity ) == 2 );
        CHECK( intersectionSize( locals.at( 1 ).connectivity, locals.at( 2 ).connectivity ) == 1 );
        CHECK( intersectionSize( locals.at( 2 ).connectivity, locals.at( 3 ).connectivity ) == 0 );

        for( size_t elem_ii = 0; elem_ii < locals.size(); elem_ii++ )
        {
            checkExtractionValues(
                kv, degree, elem_ii, locals.at( elem_ii ).connectivity, locals.at( elem_ii ).extraction, 1e-12 );
        }
    }
}

TEST_CASE( "Tensor-product element connectivity respects mixed C0 and broken continuity" )
{
    const double ptol = 1e-10;
    const KnotVector kv_s( { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 }, ptol );
    const KnotVector kv_t( { 0, 0, 1, 1, 2, 2, 3, 3 }, ptol );
    const TPSplineSpace ss = buildBSpline( { kv_s, kv_t }, { 2, 1 } );
    const auto& cmap = dynamic_cast<const TPCombinatorialMap&>( ss.basisComplex().parametricAtlas().cmap() );

    std::map<std::pair<size_t, size_t>, std::vector<FunctionId>> conn_by_elem;
    iterateCellsWhile( cmap, 2, [&]( const Cell& c ) {
        const auto [source_cell, line_cell] = unflattenCell( cmap, c );
        REQUIRE( source_cell.has_value() );
        REQUIRE( line_cell.has_value() );

        const Eigen::MatrixXd op = ss.extractionOperator( c );
        CHECK( ss.connectivity( c ).size() == 6 );
        CHECK( op.rows() == 6 );
        CHECK( op.cols() == 6 );
        CHECK( ( op.colwise().sum().transpose() - Eigen::VectorXd::Ones( op.cols() ) ).norm() < 1e-14 );

        conn_by_elem.emplace(
            std::make_pair( source_cell.value().dart().id(), line_cell.value().dart().id() ), ss.connectivity( c ) );
        return true;
    } );

    REQUIRE( conn_by_elem.size() == 9 );

    for( size_t t_elem = 0; t_elem < 3; t_elem++ )
    {
        CHECK( intersectionSize( conn_by_elem.at( { 0, t_elem } ), conn_by_elem.at( { 1, t_elem } ) ) == 2 );
        CHECK( intersectionSize( conn_by_elem.at( { 1, t_elem } ), conn_by_elem.at( { 2, t_elem } ) ) == 2 );
    }

    for( size_t s_elem = 0; s_elem < 3; s_elem++ )
    {
        CHECK( intersectionSize( conn_by_elem.at( { s_elem, 0 } ), conn_by_elem.at( { s_elem, 1 } ) ) == 0 );
        CHECK( intersectionSize( conn_by_elem.at( { s_elem, 1 } ), conn_by_elem.at( { s_elem, 2 } ) ) == 0 );
    }
}

TEST_CASE( "Local extraction helper matches BSplineSpace1d element extraction" )
{
    const double ptol = 1e-10;
    const KnotVector kv( { 0, 0, 0, 0, 1, 1, 2, 3, 3, 3, 3 }, ptol );
    const size_t degree = 3;

    const auto cmap = std::make_shared<const CombinatorialMap1d>( numElements( kv ) );
    const auto param = std::make_shared<const ParametricAtlas1d>( cmap, parametricLengths( kv ) );
    const auto bc = std::make_shared<const BasisComplex1d>( param, degree );
    const BSplineSpace1d ss( bc, kv );

    const std::vector<LocalExtraction> locals = localExtractions( kv, degree );
    REQUIRE( locals.size() == numElements( kv ) );

    iterateCellsWhile( *cmap, 1, [&]( const Cell& c ) {
        const size_t elem_ii = c.dart().id();
        const LocalExtraction& local = locals.at( elem_ii );
        const std::vector<FunctionId> ss_conn = ss.connectivity( c );

        REQUIRE( local.connectivity.size() == ss_conn.size() );
        for( size_t i = 0; i < ss_conn.size(); i++ )
        {
            CHECK( local.connectivity.at( i ) == ss_conn.at( i ) );
        }
        CHECK( ( local.extraction - ss.extractionOperator( c ) ).norm() < 1e-14 );

        checkExtractionValues( kv, degree, elem_ii, local.connectivity, local.extraction, 1e-12 );
        return true;
    } );
}

TEST_CASE( "BSplineSpace1d element extraction supports broken C-1 knot vectors" )
{
    const double ptol = 1e-10;
    const size_t degree = 2;
    const KnotVector kv( { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3 }, ptol );

    const auto cmap = std::make_shared<const CombinatorialMap1d>( numElements( kv ) );
    const auto param = std::make_shared<const ParametricAtlas1d>( cmap, parametricLengths( kv ) );
    const auto bc = std::make_shared<const BasisComplex1d>( param, degree );
    const BSplineSpace1d ss( bc, kv );

    CHECK( ss.numFunctions() == numElements( kv ) * ( degree + 1 ) );

    std::map<size_t, std::vector<FunctionId>> conn_by_elem;
    iterateCellsWhile( *cmap, 1, [&]( const Cell& c ) {
        const size_t elem_ii = c.dart().id();
        const std::vector<FunctionId> conn = ss.connectivity( c );
        const Eigen::MatrixXd op = ss.extractionOperator( c );

        CHECK( conn.size() == degree + 1 );
        CHECK( ( op - Eigen::MatrixXd::Identity( degree + 1, degree + 1 ) ).norm() < 1e-14 );
        checkExtractionValues( kv, degree, elem_ii, conn, op, 1e-12 );

        conn_by_elem.emplace( elem_ii, conn );
        return true;
    } );

    REQUIRE( conn_by_elem.size() == numElements( kv ) );
    CHECK( intersectionSize( conn_by_elem.at( 0 ), conn_by_elem.at( 1 ) ) == 0 );
    CHECK( intersectionSize( conn_by_elem.at( 1 ), conn_by_elem.at( 2 ) ) == 0 );
}
