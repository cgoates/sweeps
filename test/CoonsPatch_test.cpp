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
#include <CoonsPatch.hpp>
#include <IndexOperations.hpp>

using namespace basis;
using namespace param;
using namespace topology;

TEST_CASE( "2d coons patch square" )
{
    const KnotVector kv( {0,0,0,0,3,6,6,6,6}, 1e-10 );
    const size_t degree = 3;

    const CombinatorialMap1d cmap_1d( numElements( kv ) );
    const ParametricAtlas1d param_1d( cmap_1d, parametricLengths( kv ) );
    const BasisComplex1d bc_1d( param_1d, degree );
    const BSplineSpace1d ss_1d( bc_1d, kv );

    const TPCombinatorialMap cmap_2d( cmap_1d, cmap_1d );
    const TPParametricAtlas param_2d( cmap_2d, param_1d, param_1d );
    const TPBasisComplex bc_2d( param_2d, bc_1d, bc_1d );
    const TPSplineSpace ss_2d( bc_2d, ss_1d, ss_1d );

    const Eigen::VectorXd grev = grevillePoints( kv, degree );

    const std::vector<Eigen::MatrixXd> bdry_coeffs{
        ( Eigen::MatrixXd( grev.size(), 2 ) << Eigen::VectorXd::Zero( grev.size() ), grev ).finished(),
        ( Eigen::MatrixXd( grev.size(), 2 ) << Eigen::VectorXd::Constant( grev.size(), 6 ), grev ).finished(),
        ( Eigen::MatrixXd( grev.size(), 2 ) << grev, Eigen::VectorXd::Zero( grev.size() ) ).finished(),
        ( Eigen::MatrixXd( grev.size(), 2 ) << grev, Eigen::VectorXd::Constant( grev.size(), 6 ) ).finished()
    };

    const auto cpts = fitting::coonsPatch( ss_2d, bdry_coeffs );

    for( Eigen::Index i = 0; i < cpts.rows(); i++ )
    {
        const auto tp_id = util::unflatten( i, {ss_1d.numFunctions(), ss_1d.numFunctions() } );
        CHECK( util::equals( cpts.row( i ), Eigen::Vector2d( grev( tp_id.at( 0 ) ), grev( tp_id.at( 1 ) ) ), 1e-12 ) );
    }
}

TEST_CASE( "3d coons patch perturbed cube" )
{
    const KnotVector kv( {0,0,0,6,6,6}, 1e-10 );
    const size_t degree = 2;

    const CombinatorialMap1d cmap_1d( numElements( kv ) );
    const ParametricAtlas1d param_1d( cmap_1d, parametricLengths( kv ) );
    const BasisComplex1d bc_1d( param_1d, degree );
    const BSplineSpace1d ss_1d( bc_1d, kv );

    const TPCombinatorialMap cmap_2d( cmap_1d, cmap_1d );
    const TPParametricAtlas param_2d( cmap_2d, param_1d, param_1d );
    const TPBasisComplex bc_2d( param_2d, bc_1d, bc_1d );
    const TPSplineSpace ss_2d( bc_2d, ss_1d, ss_1d );

    const TPCombinatorialMap cmap_3d( cmap_2d, cmap_1d );
    const TPParametricAtlas param_3d( cmap_3d, param_2d, param_1d );
    const TPBasisComplex bc_3d( param_3d, bc_2d, bc_1d );
    const TPSplineSpace ss_3d( bc_3d, ss_2d, ss_1d );

    const Eigen::VectorXd grev = grevillePoints( kv, degree );

    const Eigen::MatrixX2d grev2d = util::tensorProduct( { grev, grev } );
    const Eigen::VectorXd zeros = Eigen::VectorXd::Zero( grev2d.rows() );
    const Eigen::VectorXd sixes = Eigen::VectorXd::Constant( grev2d.rows(), 6 );

    std::vector<Eigen::MatrixXd> bdry_coeffs{
        ( Eigen::MatrixXd( grev2d.rows(), 3 ) << zeros, grev2d ).finished(),
        ( Eigen::MatrixXd( grev2d.rows(), 3 ) << sixes, grev2d ).finished(),
        ( Eigen::MatrixXd( grev2d.rows(), 3 ) << grev2d.col( 0 ), zeros, grev2d.col( 1 ) ).finished(),
        ( Eigen::MatrixXd( grev2d.rows(), 3 ) << grev2d.col( 0 ), sixes, grev2d.col( 1 ) ).finished(),
        ( Eigen::MatrixXd( grev2d.rows(), 3 ) << grev2d, zeros ).finished(),
        ( Eigen::MatrixXd( grev2d.rows(), 3 ) << grev2d, sixes ).finished()
    };

    bdry_coeffs.at( 1 )( 4, 0 ) = 8;

    const auto cpts = fitting::coonsPatch( ss_3d, bdry_coeffs );

    CHECK( util::equals( cpts( 13, 0 ), 4, 1e-12 ) );
}