#include <catch2/catch_test_macros.hpp>
#include <GenericSplineSpace.hpp>
#include <BasisComplex1d.hpp>
#include <ParametricAtlas1d.hpp>
#include <CombinatorialMap1d.hpp>
#include <KnotVector.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

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