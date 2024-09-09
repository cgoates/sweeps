#include <catch2/catch_test_macros.hpp>
#include <HierarchicalTPParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

using namespace topology;
using namespace param;

TEST_CASE( "Simple hierarchical atlas 1" )
{
    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( 4 );
    const auto topo1d_8 = std::make_shared<const CombinatorialMap1d>( 8 );
    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_2 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_4 );
    const auto tp_topo_3 = std::make_shared<const TPCombinatorialMap>( topo1d_4, topo1d_8 );

    const auto atlas1d_1 = std::make_shared<const ParametricAtlas1d>( topo1d_1, Eigen::VectorXd::Ones( 1 ) );
    const auto atlas1d_2 = std::make_shared<const ParametricAtlas1d>( topo1d_2, Eigen::VectorXd::Constant( 2, 0.5 ) );
    const auto atlas1d_4 = std::make_shared<const ParametricAtlas1d>( topo1d_4, Eigen::VectorXd::Constant( 4, 0.25 ) );
    const auto atlas1d_8 = std::make_shared<const ParametricAtlas1d>( topo1d_8, Eigen::VectorXd::Constant( 8, 0.125 ) );
    const auto tp_atlas_1 = std::make_shared<const TPParametricAtlas>( tp_topo_1, atlas1d_1, atlas1d_2 );
    const auto tp_atlas_2 = std::make_shared<const TPParametricAtlas>( tp_topo_2, atlas1d_2, atlas1d_4 );
    const auto tp_atlas_3 = std::make_shared<const TPParametricAtlas>( tp_topo_3, atlas1d_4, atlas1d_8 );

    const std::shared_ptr<const HierarchicalTPCombinatorialMap> cmap( new HierarchicalTPCombinatorialMap( { tp_topo_1, tp_topo_2, tp_topo_3 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ), Face( Dart( 12 ) ) },
        {}
    } ) );

    const HierarchicalTPParametricAtlas atlas( cmap, {tp_atlas_1, tp_atlas_2, tp_atlas_3} );

    iterateCellsWhile( *cmap, 2, [&]( const Face& f ) {
        CHECK( atlas.parentDomain( f ) == cubeDomain( 2 ) );
        return true;
    } );

    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 5 ) ) ), Eigen::Vector2d( 1.0, 0.5 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 8 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 12 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 17 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 20 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 28 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 0.5, 0.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { false, false, false, true } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 24 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 0.0, 0.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { false, true, false, true } ) );
    }
}

TEST_CASE( "Simple hierarchical atlas 2" )
{
    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( 4 );
    const auto topo1d_8 = std::make_shared<const CombinatorialMap1d>( 8 );
    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_2 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_4 );
    const auto tp_topo_3 = std::make_shared<const TPCombinatorialMap>( topo1d_4, topo1d_8 );

    const auto atlas1d_1 = std::make_shared<const ParametricAtlas1d>( topo1d_1, Eigen::VectorXd::Ones( 1 ) );
    const auto atlas1d_2 = std::make_shared<const ParametricAtlas1d>( topo1d_2, Eigen::VectorXd::Constant( 2, 0.5 ) );
    const auto atlas1d_4 = std::make_shared<const ParametricAtlas1d>( topo1d_4, Eigen::VectorXd::Constant( 4, 0.25 ) );
    const auto atlas1d_8 = std::make_shared<const ParametricAtlas1d>( topo1d_8, Eigen::VectorXd::Constant( 8, 0.125 ) );
    const auto tp_atlas_1 = std::make_shared<const TPParametricAtlas>( tp_topo_1, atlas1d_1, atlas1d_2 );
    const auto tp_atlas_2 = std::make_shared<const TPParametricAtlas>( tp_topo_2, atlas1d_2, atlas1d_4 );
    const auto tp_atlas_3 = std::make_shared<const TPParametricAtlas>( tp_topo_3, atlas1d_4, atlas1d_8 );

    const std::shared_ptr<const HierarchicalTPCombinatorialMap> cmap( new HierarchicalTPCombinatorialMap( { tp_topo_1, tp_topo_2, tp_topo_3 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ) },
        { Face( Dart( 40 ) ), Face( Dart( 44 ) ), Face( Dart( 56 ) ), Face( Dart( 60 ) ) }
    } ) );

    const HierarchicalTPParametricAtlas atlas( cmap, {tp_atlas_1, tp_atlas_2, tp_atlas_3} );

    iterateCellsWhile( *cmap, 2, [&]( const Face& f ) {
        CHECK( atlas.parentDomain( f ) == cubeDomain( 2 ) );
        return true;
    } );

    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 5 ) ) ), Eigen::Vector2d( 1.0, 0.5 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 8 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 12 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 16 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 80 ) ) ), Eigen::Vector2d( 0.25, 0.125 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 84 ) ) ), Eigen::Vector2d( 0.25, 0.125 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 96 ) ) ), Eigen::Vector2d( 0.25, 0.125 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 100 ) ) ), Eigen::Vector2d( 0.25, 0.125 ), 1e-12 ) );

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 112 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 0.5, 0.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { false, false, false, true } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 116 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 0.75, 0.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { false, false, false, true } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 93 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 1.0, 0.5 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { true, false, false, false } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 66 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 0.5, 1.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { false, false, true, false } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 70 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 1.0, 1.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { true, false, true, false } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 77 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 1.0, 0.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { true, false, false, true } ) );
    }

}

TEST_CASE( "Simple hierarchical atlas 3" )
{
    const auto topo1d_1 = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto topo1d_2 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto topo1d_4 = std::make_shared<const CombinatorialMap1d>( 4 );
    const auto topo1d_8 = std::make_shared<const CombinatorialMap1d>( 8 );
    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( topo1d_1, topo1d_2 );
    const auto tp_topo_2 = std::make_shared<const TPCombinatorialMap>( topo1d_2, topo1d_4 );
    const auto tp_topo_3 = std::make_shared<const TPCombinatorialMap>( topo1d_4, topo1d_8 );

    const auto atlas1d_1 = std::make_shared<const ParametricAtlas1d>( topo1d_1, Eigen::VectorXd::Ones( 1 ) );
    const auto atlas1d_2 = std::make_shared<const ParametricAtlas1d>( topo1d_2, Eigen::VectorXd::Constant( 2, 0.5 ) );
    const auto atlas1d_4 = std::make_shared<const ParametricAtlas1d>( topo1d_4, Eigen::VectorXd::Constant( 4, 0.25 ) );
    const auto atlas1d_8 = std::make_shared<const ParametricAtlas1d>( topo1d_8, Eigen::VectorXd::Constant( 8, 0.125 ) );
    const auto tp_atlas_1 = std::make_shared<const TPParametricAtlas>( tp_topo_1, atlas1d_1, atlas1d_2 );
    const auto tp_atlas_2 = std::make_shared<const TPParametricAtlas>( tp_topo_2, atlas1d_2, atlas1d_4 );
    const auto tp_atlas_3 = std::make_shared<const TPParametricAtlas>( tp_topo_3, atlas1d_4, atlas1d_8 );

    const std::shared_ptr<const HierarchicalTPCombinatorialMap> cmap( new HierarchicalTPCombinatorialMap( { tp_topo_1, tp_topo_2, tp_topo_3 }, {
        { Face( Dart( 0 ) ) },
        { Face( Dart( 16 ) ), Face( Dart( 24 ) ), Face( Dart( 28 ) ) },
        { Face( Dart( 72 ) ), Face( Dart( 76 ) ), Face( Dart( 88 ) ), Face( Dart( 92 ) ) }
    } ) );

    const HierarchicalTPParametricAtlas atlas( cmap, {tp_atlas_1, tp_atlas_2, tp_atlas_3} );

    iterateCellsWhile( *cmap, 2, [&]( const Face& f ) {
        CHECK( atlas.parentDomain( f ) == cubeDomain( 2 ) );
        return true;
    } );

    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 0 ) ) ), Eigen::Vector2d( 1.0, 0.5 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 24 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 32 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 112 ) ) ), Eigen::Vector2d( 0.25, 0.125 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 116 ) ) ), Eigen::Vector2d( 0.25, 0.125 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 128 ) ) ), Eigen::Vector2d( 0.25, 0.125 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 132 ) ) ), Eigen::Vector2d( 0.25, 0.125 ), 1e-12 ) );
    CHECK( util::equals( atlas.parametricLengths( Face( Dart( 144 ) ) ), Eigen::Vector2d( 0.5, 0.25 ), 1e-12 ) );

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 102 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 1.0, 1.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { true, false, true, false } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 98 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 0.75, 1.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { false, false, true, false } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 18 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 0.5, 1.0 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { false, false, true, false } ) );
    }

    {
        const ParentPoint ppt = atlas.parentPoint( Vertex( Dart( 125 ) ) );
        CHECK( ppt.mDomain == cubeDomain( 2 ) );
        CHECK( util::equals( ppt.mPoint, Eigen::Vector2d( 1.0, 0.5 ), 1e-9 ) );
        CHECK( ppt.mBaryCoordIsZero == BaryCoordIsZeroVec( { true, false, false, false } ) );
    }

}