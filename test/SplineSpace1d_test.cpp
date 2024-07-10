#include <catch2/catch_test_macros.hpp>
#include <GenericSplineSpace.hpp>
#include <BasisComplex1d.hpp>
#include <ParametricAtlas1d.hpp>
#include <CombinatorialMap1d.hpp>
#include <KnotVector.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <TPSplineSpace.hpp>

using namespace basis;
using namespace param;
using namespace topology;

TEST_CASE( "Extraction of cubic spline with various knot multiplicities" )
{
    const KnotVector kv( {0,0,0,0,1,1,2,3,3,3,3}, 1e-10 );
    const size_t degree = 3;

    const CombinatorialMap1d cmap( 3 );
    const ParametricAtlas1d param( cmap, parametricLengths( kv ) );
    const BasisComplex1d bc( param, degree );
    const GenericSplineSpace ss = knotVectorSplineSpace( bc, kv );

    iterateCellsWhile( cmap, 1, [&]( const topology::Edge& e ) {
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

    const CombinatorialMap1d cmap( 2 );
    const ParametricAtlas1d param( cmap, parametricLengths( kv1 ) );
    const BasisComplex1d bc1( param, degree1 );
    const BasisComplex1d bc2( param, degree2 );
    const GenericSplineSpace ss1 = knotVectorSplineSpace( bc1, kv1 );
    const GenericSplineSpace ss2 = knotVectorSplineSpace( bc2, kv2 );

    const TPCombinatorialMap cmap_2d( cmap, cmap );
    const TPParametricAtlas param_2d( cmap_2d, param, param );
    const TPBasisComplex bc_2d( param_2d, bc1, bc2 );
    const TPSplineSpace ss_2d( bc_2d, ss1, ss2 );

    iterateCellsWhile( cmap_2d, 2, [&]( const topology::Face& f ) {
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