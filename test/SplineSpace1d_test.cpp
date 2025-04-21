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
#include <Eigen/Dense>

using namespace basis;
using namespace param;
using namespace topology;

TEST_CASE( "Extraction of cubic spline with various knot multiplicities" )
{
    const KnotVector kv( {0,0,0,0,1,1,2,3,3,3,3}, 1e-10 );
    const size_t degree = 3;

    const auto cmap = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto param = std::make_shared<const ParametricAtlas1d>( cmap, parametricLengths( kv ) );
    const auto bc = std::make_shared<const BasisComplex1d>( param, degree );
    const BSplineSpace1d ss( bc, kv );

    iterateCellsWhile( *cmap, 1, [&]( const topology::Edge& e ) {
        const Eigen::MatrixXd op = ss.extractionOperator( e );
        CHECK( util::equals( op.colwise().sum(), Eigen::VectorXd::Ones( op.cols() ), 1e-12 ) );
        CHECK( ss.connectivity( e ).size() == degree + 1 );
        CHECK( op.rows() == degree + 1 );
        CHECK( op.rows() == op.cols() );

        Eigen::FullPivLU<Eigen::MatrixXd> lu( op );
        CHECK( lu.rank() == op.rows() );
        return true;
    } );
}

TEST_CASE( "TP Spline space" )
{
    const KnotVector kv1( {0,0,0,0,1,2,2,2,2}, 1e-10 );
    const KnotVector kv2( {0,0,0,1,2,2,2}, 1e-10 );
    const size_t degree1 = 3;
    const size_t degree2 = 2;

    const auto cmap = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto param = std::make_shared<const ParametricAtlas1d>( cmap, parametricLengths( kv1 ) );
    const auto bc1 = std::make_shared<const BasisComplex1d>( param, degree1 );
    const auto bc2 = std::make_shared<const BasisComplex1d>( param, degree2 );
    const auto ss1 = std::make_shared<const BSplineSpace1d>( bc1, kv1 );
    const auto ss2 = std::make_shared<const BSplineSpace1d>( bc2, kv2 );

    const auto cmap_2d = std::make_shared<const TPCombinatorialMap>( cmap, cmap );
    const auto param_2d = std::make_shared<const TPParametricAtlas>( cmap_2d, param, param );
    const auto bc_2d = std::make_shared<const TPBasisComplex>( param_2d, bc1, bc2 );
    const TPSplineSpace ss_2d( bc_2d, ss1, ss2 );

    iterateCellsWhile( *cmap_2d, 2, [&]( const topology::Face& f ) {
        const Eigen::MatrixXd op = ss_2d.extractionOperator( f );
        CHECK( util::equals( op.colwise().sum(), Eigen::VectorXd::Ones( op.cols() ), 1e-12 ) );
        CHECK( ss_2d.connectivity( f ).size() == ( degree1 + 1 ) * ( degree2 + 1 ) );
        CHECK( op.rows() == ( degree1 + 1 ) * ( degree2 + 1 ) );
        CHECK( op.rows() == op.cols() );

        Eigen::FullPivLU<Eigen::MatrixXd> lu( op );
        CHECK( lu.rank() == op.rows() );
        return true;
    } );
}