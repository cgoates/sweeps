#include <catch2/catch_test_macros.hpp>
#include <HierarchicalTPSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <VTKOutput.hpp>
#include <SplineSpaceEvaluator.hpp>

using namespace topology;
using namespace param;
using namespace basis;

constexpr bool VTK_OUTPUT = false;

TEST_CASE( "Simple hierarchical atlas 1" )
{
    const KnotVector kv1( {0, 0, 0, 1, 1, 1}, 1e-10 );
    const KnotVector kv2( {0, 0, 0, 0.5, 1, 1, 1}, 1e-10 );
    const KnotVector kv4( {0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1}, 1e-10 );

    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( numElements( kv1 ) );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( numElements( kv2 ) );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( numElements( kv4 ) );
    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_2 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_4 );

    const auto atlas1d_1 = std::make_shared<const ParametricAtlas1d>( topo1d_1, parametricLengths( kv1 ) );
    const auto atlas1d_2 = std::make_shared<const ParametricAtlas1d>( topo1d_2, parametricLengths( kv2 ) );
    const auto atlas1d_4 = std::make_shared<const ParametricAtlas1d>( topo1d_4, parametricLengths( kv4 ) );
    const auto tp_atlas_1 = std::make_shared<const TPParametricAtlas>( tp_topo_1, atlas1d_1, atlas1d_2 );
    const auto tp_atlas_2 = std::make_shared<const TPParametricAtlas>( tp_topo_2, atlas1d_2, atlas1d_4 );

    const size_t degree = 2;
    const auto bc1d_1 = std::make_shared<const BasisComplex1d>( atlas1d_1, degree );
    const auto bc1d_2 = std::make_shared<const BasisComplex1d>( atlas1d_2, degree );
    const auto bc1d_4 = std::make_shared<const BasisComplex1d>( atlas1d_4, degree );
    const auto tp_bc_1 = std::make_shared<const TPBasisComplex>( tp_atlas_1, bc1d_1, bc1d_2 );
    const auto tp_bc_2 = std::make_shared<const TPBasisComplex>( tp_atlas_2, bc1d_2, bc1d_4 );

    const auto ss1d_1 = std::make_shared<const BSplineSpace1d>( bc1d_1, kv1 );
    const auto ss1d_2 = std::make_shared<const BSplineSpace1d>( bc1d_2, kv2 );
    const auto ss1d_4 = std::make_shared<const BSplineSpace1d>( bc1d_4, kv4 );
    const auto tp_ss_1 = std::make_shared<const TPSplineSpace>( tp_bc_1, ss1d_1, ss1d_2 );
    const auto tp_ss_2 = std::make_shared<const TPSplineSpace>( tp_bc_2, ss1d_2, ss1d_4 );

    const std::shared_ptr<const HierarchicalTPCombinatorialMap> cmap( new HierarchicalTPCombinatorialMap( { tp_topo_1, tp_topo_2 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ), Face( Dart( 12 ) ) }
    } ) );

    const std::shared_ptr<const HierarchicalTPParametricAtlas> atlas(
        new HierarchicalTPParametricAtlas( cmap, { tp_atlas_1, tp_atlas_2 } ) );

    const std::shared_ptr<const HierarchicalTPBasisComplex> bc(
        new HierarchicalTPBasisComplex( atlas, { tp_bc_1, tp_bc_2 } ) );

    const HierarchicalTPSplineSpace ss( bc, { tp_ss_1, tp_ss_2 }, { {3, 4, 5, 6, 7, 8, 9, 10, 11}, {0, 1, 2, 3, 4, 5, 6, 7} } );

    const Eigen::MatrixX3d geom = ( Eigen::MatrixX3d( 17, 3 ) <<
        0, 0.5, 0,
        0.5, 0.5, 0,
        1, 0.5, 0,
        0, 1.5, 0,
        0.5, 1.5, 0,
        1, 1.5, 0,
        0, 2, 0,
        0.5, 2, 0,
        1, 2, 0,
        0, 0, 0,
        0.25, 0, 0,
        0.75, 0, 0,
        1, 0, 0,
        0, 0.25, 0,
        0.25, 0.25, 0,
        0.75, 0.25, 0,
        1, 0.25, 0 ).finished();

    if( VTK_OUTPUT ) io::outputBezierMeshToVTK( ss, geom, "hier_bez_test.vtu" );

    eval::SplineSpaceEvaluator evaler( ss, 0 );

    const auto test_center_pt = [&]( const Face& f, const Eigen::Vector3d& expected ) {
        const ParentPoint ppt = pointOnBoundary( atlas->parentDomain( f ), parentDomainBoundary( *atlas, f ) );
        evaler.localizeElement( f );
        evaler.localizePoint( ppt );

        const Eigen::Vector3d pt = evaler.evaluateManifold( geom.transpose() );

        CHECK( util::equals( pt, expected, 1e-9 ) );
    };

    test_center_pt( Face( 5 ), Eigen::Vector3d( 0.5, 1.5, 0.0 ) );
    test_center_pt( Face( 8 ), Eigen::Vector3d( 0.25, 0.25, 0.0 ) );
    test_center_pt( Face( 12 ), Eigen::Vector3d( 0.75, 0.25, 0.0 ) );
}

TEST_CASE( "Simple hierarchical atlas 2" )
{
    const KnotVector kv1( {0, 0, 1, 1}, 1e-10 );
    const KnotVector kv2( {0, 0, 0.5, 1, 1}, 1e-10 );
    const KnotVector kv4( {0, 0, 0.25, 0.5, 0.75, 1, 1}, 1e-10 );
    const KnotVector kv8( { 0, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1 }, 1e-10 );

    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( numElements( kv1 ) );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( numElements( kv2 ) );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( numElements( kv4 ) );
    const auto topo1d_8 = std::make_shared<const CombinatorialMap1d>( numElements( kv8 ) );
    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_2 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_4 );
    const auto tp_topo_3 = std::make_shared<const TPCombinatorialMap>( topo1d_4, topo1d_8 );

    const auto atlas1d_1 = std::make_shared<const ParametricAtlas1d>( topo1d_1, parametricLengths( kv1 ) );
    const auto atlas1d_2 = std::make_shared<const ParametricAtlas1d>( topo1d_2, parametricLengths( kv2 ) );
    const auto atlas1d_4 = std::make_shared<const ParametricAtlas1d>( topo1d_4, parametricLengths( kv4 ) );
    const auto atlas1d_8 = std::make_shared<const ParametricAtlas1d>( topo1d_8, parametricLengths( kv8 ) );
    const auto tp_atlas_1 = std::make_shared<const TPParametricAtlas>( tp_topo_1, atlas1d_1, atlas1d_2 );
    const auto tp_atlas_2 = std::make_shared<const TPParametricAtlas>( tp_topo_2, atlas1d_2, atlas1d_4 );
    const auto tp_atlas_3 = std::make_shared<const TPParametricAtlas>( tp_topo_3, atlas1d_4, atlas1d_8 );

    const size_t degree = 1;
    const auto bc1d_1 = std::make_shared<const BasisComplex1d>( atlas1d_1, degree );
    const auto bc1d_2 = std::make_shared<const BasisComplex1d>( atlas1d_2, degree );
    const auto bc1d_4 = std::make_shared<const BasisComplex1d>( atlas1d_4, degree );
    const auto bc1d_8 = std::make_shared<const BasisComplex1d>( atlas1d_8, degree );
    const auto tp_bc_1 = std::make_shared<const TPBasisComplex>( tp_atlas_1, bc1d_1, bc1d_2 );
    const auto tp_bc_2 = std::make_shared<const TPBasisComplex>( tp_atlas_2, bc1d_2, bc1d_4 );
    const auto tp_bc_3 = std::make_shared<const TPBasisComplex>( tp_atlas_3, bc1d_4, bc1d_8 );

    const auto ss1d_1 = std::make_shared<const BSplineSpace1d>( bc1d_1, kv1 );
    const auto ss1d_2 = std::make_shared<const BSplineSpace1d>( bc1d_2, kv2 );
    const auto ss1d_4 = std::make_shared<const BSplineSpace1d>( bc1d_4, kv4 );
    const auto ss1d_8 = std::make_shared<const BSplineSpace1d>( bc1d_8, kv8 );
    const auto tp_ss_1 = std::make_shared<const TPSplineSpace>( tp_bc_1, ss1d_1, ss1d_2 );
    const auto tp_ss_2 = std::make_shared<const TPSplineSpace>( tp_bc_2, ss1d_2, ss1d_4 );
    const auto tp_ss_3 = std::make_shared<const TPSplineSpace>( tp_bc_3, ss1d_4, ss1d_8 );

    const std::shared_ptr<const HierarchicalTPCombinatorialMap> cmap( new HierarchicalTPCombinatorialMap( { tp_topo_1, tp_topo_2, tp_topo_3 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ) },
        { Face( Dart( 40 ) ), Face( Dart( 44 ) ), Face( Dart( 56 ) ), Face( Dart( 60 ) ) }
    } ) );

    const std::shared_ptr<const HierarchicalTPParametricAtlas> atlas(
        new HierarchicalTPParametricAtlas( cmap, { tp_atlas_1, tp_atlas_2, tp_atlas_3 } ) );

    const std::shared_ptr<const HierarchicalTPBasisComplex> bc(
        new HierarchicalTPBasisComplex( atlas, { tp_bc_1, tp_bc_2, tp_bc_3 } ) );

    const HierarchicalTPSplineSpace ss( bc, { tp_ss_1, tp_ss_2, tp_ss_3 }, { {2, 3, 4, 5}, {0, 1, 2, 3, 4, 5}, {18, 19} } );

    const Eigen::MatrixX3d geom = ( Eigen::MatrixX3d( 12, 3 ) <<
        0, 1, 0,
        1, 1, 0,
        0, 2, 0,
        1, 2, 0,
        0, 0, 0,
        0.5, 0, 0,
        1, 0, 0,
        0, 0.5, 0,
        0.5, 0.5, 0,
        1, 0.5, 0,
        0.75, 0.75, 0,
        1.0, 0.75, 0 ).finished();

    if( VTK_OUTPUT ) io::outputBezierMeshToVTK( ss, geom, "hier_bez_test.vtu" );

    eval::SplineSpaceEvaluator evaler( ss, 0 );

    const auto test_center_pt = [&]( const Face& f, const Eigen::Vector3d& expected ) {
        const ParentPoint ppt = pointOnBoundary( atlas->parentDomain( f ), parentDomainBoundary( *atlas, f ) );
        evaler.localizeElement( f );
        evaler.localizePoint( ppt );

        const Eigen::Vector3d pt = evaler.evaluateManifold( geom.transpose() );

        CHECK( util::equals( pt, expected, 1e-9 ) );
    };

    test_center_pt( Face( 5 ), Eigen::Vector3d( 0.5, 1.5, 0.0 ) );
    test_center_pt( Face( 8 ), Eigen::Vector3d( 0.25, 0.25, 0.0 ) );
    test_center_pt( Face( 12 ), Eigen::Vector3d( 0.75, 0.25, 0.0 ) );
    test_center_pt( Face( 80 ), Eigen::Vector3d( 0.625, 0.625, 0.0 ) );
}