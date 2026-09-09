#include <catch2/catch_test_macros.hpp>
#include <BSplineSpace1d.hpp>
#include <BasisComplex1d.hpp>
#include <ParametricAtlas1d.hpp>
#include <CombinatorialMap1d.hpp>
#include <KnotVector.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CommonUtils.hpp>
#include <TPSplineSpace.hpp>
#include <VectorConformingBasisComplex.hpp>
#include <VectorConformingTPSplineSpace.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <set>

using namespace basis;
using namespace param;
using namespace topology;

// Phase 2: single-patch reduced-continuity de Rham spline spaces.
//
// A degree-p primal knot with interior multiplicity p is C0 (value continuous,
// first derivative discontinuous). In the div/curl-conforming complex built from
// it, the NORMAL trace of H(div) (tangential trace of H(curl)) stays continuous
// across that knot while the complementary trace drops to C-1, and the terminal
// L2 space becomes C-1 (fully discontinuous). These tests assert that continuity
// structure directly -- no weak-form jump terms, no convergence study.
//
// Geometry is the identity map (Greville control points), so det J = 1 and the
// spatial push-forwards reduce to the parametric components; the continuity
// structure is what is under test, not the Piola scaling (covered elsewhere).

namespace
{
    constexpr double TOL = 1e-11;

    ParentPoint parentPoint2d( const Eigen::Vector2d& coordinates )
    {
        return ParentPoint( cubeDomain( 2 ), coordinates,
                            { false, false, false, false } );
    }

    // The two faces of a 2-elements-in-s, 1-element-in-t patch, ordered by
    // ascending parametric s-start (left element first).
    std::array<Face, 2> sortedInterfaceElements( const TPSplineSpace& h1 )
    {
        std::vector<Face> elements;
        iterateCellsWhile(
            h1.basisComplex().parametricAtlas().cmap(), 2,
            [&]( const Face& f ) { elements.push_back( f ); return true; } );
        REQUIRE( elements.size() == 2 );
        std::sort( elements.begin(), elements.end(),
            [&]( const Face& a, const Face& b ) {
                return h1.basisComplex().parametricAtlas().parametricStarts( a )( 0 ) <
                       h1.basisComplex().parametricAtlas().parametricStarts( b )( 0 );
            } );
        return { elements.at( 0 ), elements.at( 1 ) };
    }

    // Identity-map control points (2 x numFunctions) laid out to match
    // component ordering column = j * (#s functions) + i.
    Eigen::MatrixXd grevilleIdentityGeometry(
        const KnotVector& kv_s, const size_t p_s,
        const KnotVector& kv_t, const size_t p_t,
        const size_t num_functions )
    {
        const Eigen::VectorXd gs = grevillePoints( kv_s, p_s );
        const Eigen::VectorXd gt = grevillePoints( kv_t, p_t );
        Eigen::MatrixXd cp( 2, num_functions );
        for( Eigen::Index j = 0; j < gt.size(); ++j )
            for( Eigen::Index i = 0; i < gs.size(); ++i )
            {
                const Eigen::Index col = j * gs.size() + i;
                cp( 0, col ) = gs( i );
                cp( 1, col ) = gt( j );
            }
        return cp;
    }

    // Connectivity ids whose component `component` is nonzero at the evaluated point.
    std::set<size_t> nonzeroSupport(
        const std::vector<FunctionId>& connectivity,
        const Eigen::MatrixXd& values, const Eigen::Index component )
    {
        std::set<size_t> support;
        for( size_t f = 0; f < connectivity.size(); ++f )
            if( std::abs( values( f, component ) ) > TOL )
                support.insert( connectivity.at( f ) );
        return support;
    }
}

TEST_CASE( "Repeated interior knot yields a C0 H1 space" )
{
    const KnotVector kv_t( { 0, 0, 0, 2, 2, 2 }, 1e-10 );
    const KnotVector kv_s_smooth( { 0, 0, 0, 1, 2, 2, 2 }, 1e-10 );    // C1 at s=1
    const KnotVector kv_s_reduced( { 0, 0, 0, 1, 1, 2, 2, 2 }, 1e-10 ); // C0 at s=1

    const TPSplineSpace smooth = buildBSpline( { kv_s_smooth, kv_t }, { 2, 2 } );
    const TPSplineSpace reduced = buildBSpline( { kv_s_reduced, kv_t }, { 2, 2 } );

    // Repeating the interior knot adds exactly one s-function per t-function.
    CHECK( reduced.numFunctions() == smooth.numFunctions() + 3 );

    const Eigen::Vector2d left_parent( 1.0, 0.37 );  // left element, right edge
    const Eigen::Vector2d right_parent( 0.0, 0.37 ); // right element, left edge

    const auto evaluate = [&]( const TPSplineSpace& space ) {
        eval::SplineSpaceEvaluator evals( space, 1 );
        const auto sides = sortedInterfaceElements( space );

        evals.localizeElement( sides.at( 0 ) );
        evals.localizeParentPoint( parentPoint2d( left_parent ) );
        const std::vector<FunctionId> left_conn = space.connectivity( sides.at( 0 ) );
        const Eigen::MatrixXd left_val = evals.evaluateBasisValuesAtParentPoint();
        const Eigen::MatrixXd left_der =
            evals.evaluateBasisFirstDerivativesWrtParametricCoordinates();

        evals.localizeElement( sides.at( 1 ) );
        evals.localizeParentPoint( parentPoint2d( right_parent ) );
        const std::vector<FunctionId> right_conn = space.connectivity( sides.at( 1 ) );
        const Eigen::MatrixXd right_val = evals.evaluateBasisValuesAtParentPoint();
        const Eigen::MatrixXd right_der =
            evals.evaluateBasisFirstDerivativesWrtParametricCoordinates();

        // For every function shared across the interface, compare value and the
        // s-derivative (column 0 of the parametric-frame first derivatives).
        size_t shared = 0;
        size_t value_jumps = 0;
        size_t derivative_jumps = 0;
        for( size_t l = 0; l < left_conn.size(); ++l )
        {
            const auto it = std::find( right_conn.begin(), right_conn.end(), left_conn.at( l ) );
            if( it == right_conn.end() ) continue;
            const size_t r = std::distance( right_conn.begin(), it );
            ++shared;
            if( std::abs( left_val( l, 0 ) - right_val( r, 0 ) ) > TOL ) ++value_jumps;
            if( std::abs( left_der( l, 0 ) - right_der( r, 0 ) ) > TOL ) ++derivative_jumps;
        }
        return std::make_tuple( shared, value_jumps, derivative_jumps );
    };

    const auto [smooth_shared, smooth_value_jumps, smooth_derivative_jumps] = evaluate( smooth );
    const auto [reduced_shared, reduced_value_jumps, reduced_derivative_jumps] = evaluate( reduced );

    // Smooth (C1): value and first derivative both continuous across the knot.
    CHECK( smooth_shared > 0 );
    CHECK( smooth_value_jumps == 0 );
    CHECK( smooth_derivative_jumps == 0 );

    // Reduced (C0): value still continuous, but the first derivative jumps.
    CHECK( reduced_shared > 0 );
    CHECK( reduced_value_jumps == 0 );
    CHECK( reduced_derivative_jumps > 0 );
}

TEST_CASE( "Reduced-continuity H(div) preserves the normal trace, breaks the tangential" )
{
    const KnotVector kv_s( { 0, 0, 0, 1, 1, 2, 2, 2 }, 1e-10 ); // C0 at s=1
    const KnotVector kv_t( { 0, 0, 0, 2, 2, 2 }, 1e-10 );
    const TPSplineSpace h1 = buildBSpline( { kv_s, kv_t }, { 2, 2 } );

    const auto hdiv_complex = std::make_shared<const VectorConformingBasisComplex>(
        h1.basisComplexPtr(), ConformingType::Divergence );
    const VectorConformingTPSplineSpace hdiv( hdiv_complex, h1 );

    const Eigen::MatrixXd cp =
        grevilleIdentityGeometry( kv_s, 2, kv_t, 2, h1.numFunctions() );

    const auto sides = sortedInterfaceElements( h1 );
    eval::SplineSpaceEvaluator geom_evals( h1, 1 );
    eval::SplineSpaceEvaluator hdiv_evals( hdiv, 0 );

    const auto evaluate_side = [&]( const Face& elem, const Eigen::Vector2d& pp ) {
        geom_evals.localizeElement( elem );
        hdiv_evals.localizeElement( elem );
        geom_evals.localizeParentPoint( parentPoint2d( pp ) );
        hdiv_evals.localizeParentPoint( parentPoint2d( pp ) );
        return std::make_pair(
            hdiv.connectivity( elem ),
            eval::evaluateSpatialHDivBasisValues( hdiv_evals, geom_evals, cp ) );
    };

    // Check along the whole interface, not a single point.
    for( const double t : { 0.2, 0.5, 0.8 } )
    {
        const auto [left_conn, left_val] = evaluate_side( sides.at( 0 ), { 1.0, t } );
        const auto [right_conn, right_val] = evaluate_side( sides.at( 1 ), { 0.0, t } );

        // Normal component (index 0, the s-component) is C0: shared functions
        // match, and at least one carries a nonzero normal trace.
        size_t matched_normal = 0;
        for( size_t l = 0; l < left_conn.size(); ++l )
        {
            const auto it = std::find( right_conn.begin(), right_conn.end(), left_conn.at( l ) );
            if( it == right_conn.end() ) continue;
            const size_t r = std::distance( right_conn.begin(), it );
            CHECK( std::abs( left_val( l, 0 ) - right_val( r, 0 ) ) < TOL );
            if( std::abs( left_val( l, 0 ) ) > TOL ) ++matched_normal;
        }
        CHECK( matched_normal > 0 );

        // Tangential component (index 1, the t-component) is C-1: its support does
        // not bridge the interface -- no id is tangentially active on both sides.
        const std::set<size_t> left_tangential = nonzeroSupport( left_conn, left_val, 1 );
        const std::set<size_t> right_tangential = nonzeroSupport( right_conn, right_val, 1 );
        CHECK( not left_tangential.empty() );
        CHECK( not right_tangential.empty() );
        std::vector<size_t> shared_tangential;
        std::set_intersection(
            left_tangential.begin(), left_tangential.end(),
            right_tangential.begin(), right_tangential.end(),
            std::back_inserter( shared_tangential ) );
        CHECK( shared_tangential.empty() );
    }
}

TEST_CASE( "Reduced-continuity H(curl) preserves the tangential trace, breaks the normal" )
{
    const KnotVector kv_s( { 0, 0, 0, 1, 1, 2, 2, 2 }, 1e-10 ); // C0 at s=1
    const KnotVector kv_t( { 0, 0, 0, 2, 2, 2 }, 1e-10 );
    const TPSplineSpace h1 = buildBSpline( { kv_s, kv_t }, { 2, 2 } );

    const auto hcurl_complex = std::make_shared<const VectorConformingBasisComplex>(
        h1.basisComplexPtr(), ConformingType::Curl );
    const VectorConformingTPSplineSpace hcurl( hcurl_complex, h1 );

    const Eigen::MatrixXd cp =
        grevilleIdentityGeometry( kv_s, 2, kv_t, 2, h1.numFunctions() );

    const auto sides = sortedInterfaceElements( h1 );
    eval::SplineSpaceEvaluator geom_evals( h1, 1 );
    eval::SplineSpaceEvaluator hcurl_evals( hcurl, 0 );

    const auto evaluate_side = [&]( const Face& elem, const Eigen::Vector2d& pp ) {
        geom_evals.localizeElement( elem );
        hcurl_evals.localizeElement( elem );
        geom_evals.localizeParentPoint( parentPoint2d( pp ) );
        hcurl_evals.localizeParentPoint( parentPoint2d( pp ) );
        return std::make_pair(
            hcurl.connectivity( elem ),
            eval::evaluateSpatialHCurlBasisValues( hcurl_evals, geom_evals, cp ) );
    };

    for( const double t : { 0.2, 0.5, 0.8 } )
    {
        const auto [left_conn, left_val] = evaluate_side( sides.at( 0 ), { 1.0, t } );
        const auto [right_conn, right_val] = evaluate_side( sides.at( 1 ), { 0.0, t } );

        // Tangential component (index 1, the t-component) is C0: continuous.
        size_t matched_tangential = 0;
        for( size_t l = 0; l < left_conn.size(); ++l )
        {
            const auto it = std::find( right_conn.begin(), right_conn.end(), left_conn.at( l ) );
            if( it == right_conn.end() ) continue;
            const size_t r = std::distance( right_conn.begin(), it );
            CHECK( std::abs( left_val( l, 1 ) - right_val( r, 1 ) ) < TOL );
            if( std::abs( left_val( l, 1 ) ) > TOL ) ++matched_tangential;
        }
        CHECK( matched_tangential > 0 );

        // Normal component (index 0, the s-component) is C-1: support does not bridge.
        const std::set<size_t> left_normal = nonzeroSupport( left_conn, left_val, 0 );
        const std::set<size_t> right_normal = nonzeroSupport( right_conn, right_val, 0 );
        CHECK( not left_normal.empty() );
        CHECK( not right_normal.empty() );
        std::vector<size_t> shared_normal;
        std::set_intersection(
            left_normal.begin(), left_normal.end(),
            right_normal.begin(), right_normal.end(),
            std::back_inserter( shared_normal ) );
        CHECK( shared_normal.empty() );
    }
}

TEST_CASE( "A C0 primal yields a C-1 terminal L2 space" )
{
    const KnotVector kv_t( { 0, 0, 0, 2, 2, 2 }, 1e-10 );
    const KnotVector kv_s_smooth( { 0, 0, 0, 1, 2, 2, 2 }, 1e-10 );    // C1 primal -> C0 L2
    const KnotVector kv_s_reduced( { 0, 0, 0, 1, 1, 2, 2, 2 }, 1e-10 ); // C0 primal -> C-1 L2

    // Build the terminal L2 space (degree reduced in every direction) the same
    // way the de Rham complex does, then read its continuity at the s=1 knot.
    const auto build_l2_support = [&]( const KnotVector& kv_s ) {
        const TPSplineSpace h1 = buildBSpline( { kv_s, kv_t }, { 2, 2 } );
        const auto hdiv_complex = std::make_shared<const VectorConformingBasisComplex>(
            h1.basisComplexPtr(), ConformingType::Divergence );
        const VectorConformingTPSplineSpace hdiv( hdiv_complex, h1 );
        const auto l2_complex = std::make_shared<const TPBasisComplex>(
            h1.basisComplex().parametricAtlasPtr(),
            hdiv.reducedDegree1dBases().at( 0 )->basisComplexPtr(),
            hdiv.reducedDegree1dBases().at( 1 )->basisComplexPtr() );
        const TPSplineSpace l2( l2_complex,
            hdiv.reducedDegree1dBases().at( 0 ),
            hdiv.reducedDegree1dBases().at( 1 ) );

        const Eigen::MatrixXd cp =
            grevilleIdentityGeometry( kv_s, 2, kv_t, 2, h1.numFunctions() );

        const auto sides = sortedInterfaceElements( h1 );
        eval::SplineSpaceEvaluator geom_evals( h1, 1 );
        eval::SplineSpaceEvaluator l2_evals( l2, 0 );

        const auto side = [&]( const Face& elem, const Eigen::Vector2d& pp ) {
            geom_evals.localizeElement( elem );
            l2_evals.localizeElement( elem );
            geom_evals.localizeParentPoint( parentPoint2d( pp ) );
            l2_evals.localizeParentPoint( parentPoint2d( pp ) );
            return std::make_pair(
                l2.connectivity( elem ),
                eval::evaluateSpatialL2BasisValues( l2_evals, geom_evals, cp ) );
        };

        const auto [left_conn, left_val] = side( sides.at( 0 ), { 1.0, 0.37 } );
        const auto [right_conn, right_val] = side( sides.at( 1 ), { 0.0, 0.37 } );

        // Count ids that are active (nonzero) on both sides of the interface.
        const std::set<size_t> left_support = nonzeroSupport( left_conn, left_val, 0 );
        const std::set<size_t> right_support = nonzeroSupport( right_conn, right_val, 0 );
        std::vector<size_t> bridging;
        std::set_intersection(
            left_support.begin(), left_support.end(),
            right_support.begin(), right_support.end(),
            std::back_inserter( bridging ) );
        return bridging.size();
    };

    // C1 primal: terminal L2 is C0 -- at least one basis function bridges the knot.
    CHECK( build_l2_support( kv_s_smooth ) > 0 );
    // C0 primal: terminal L2 is C-1 -- no basis function bridges the knot.
    CHECK( build_l2_support( kv_s_reduced ) == 0 );
}
