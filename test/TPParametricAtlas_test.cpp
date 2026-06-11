#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <TPParametricAtlas.hpp>
#include <MultiPatchParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

using namespace topology;
using namespace param;

namespace
{
    Cell topCellAt( const TPCombinatorialMap& cmap, const util::IndexVec& element_indices )
    {
        FullyUnflattenedDart unflat;
        for( const size_t element_index : element_indices )
        {
            unflat.unflat_darts.push_back( Dart( element_index ) );
        }
        unflat.dart_pos =
            SmallVector<TPCombinatorialMap::TPDartPos, 2>(
                element_indices.size() - 1, TPCombinatorialMap::TPDartPos::DartPos0 );
        return Cell( flattenFull( cmap, unflat ), cmap.dim() );
    }
}

TEST_CASE( "3x4 TP combinatorial map" )
{
    const auto source_topo = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto line_topo = std::make_shared<const CombinatorialMap1d>( 4 );
    const auto tp_topo = std::make_shared<const TPCombinatorialMap>( source_topo, line_topo );

    const auto source = std::make_shared<const ParametricAtlas1d>( source_topo, Eigen::Vector3d( 1, 2, 3 ) );
    const auto line = std::make_shared<const ParametricAtlas1d>( line_topo, Eigen::Vector4d( 0.1, 0.2, 0.3, 0.4 ) );
    const TPParametricAtlas tp( tp_topo, source, line );

    const ParentDomain expected_domain = cubeDomain( 2 );
    iterateCellsWhile( *tp_topo, 2, [&]( const Cell& c ) {
        CHECK( tp.parentDomain( c ) == expected_domain );
        return true;
    } );

    std::array<ParentPoint, 4> expected_points( {
        ParentPoint( expected_domain, Eigen::Vector2d( 0.0, 0.0 ), { false, true, false, true } ),
        ParentPoint( expected_domain, Eigen::Vector2d( 1.0, 0.0 ), { true, false, false, true } ),
        ParentPoint( expected_domain, Eigen::Vector2d( 1.0, 1.0 ), { true, false, true, false } ),
        ParentPoint( expected_domain, Eigen::Vector2d( 0.0, 1.0 ), { false, true, true, false } )
    } );
    iterateDartsWhile( *tp_topo, [&]( const Dart& d ) {
        const auto ppt = tp.parentPoint( d );
        const auto& expected_ppt = expected_points.at( d.id() % 4 );
        CHECK( ppt.mDomain == expected_domain );
        CHECK( ppt.mBaryCoordIsZero == expected_ppt.mBaryCoordIsZero );
        CHECK( util::equals( ppt.mPoint, expected_ppt.mPoint, 1e-10 ) );
        return true;
    } );

    CHECK( util::equals( tp.parametricLengths( Cell( Dart( 0 ), 2 ) ), Eigen::Vector2d( 1.0, 0.1 ), 1e-10 ) );
    CHECK( util::equals( tp.parametricLengths( Cell( Dart( 20 ), 2 ) ), Eigen::Vector2d( 3.0, 0.2 ), 1e-10 ) );
    CHECK( util::equals( tp.parametricLengths( Cell( Dart( 37 ), 2 ) ), Eigen::Vector2d( 1.0, 0.4 ), 1e-10 ) );
    CHECK( util::equals( tp.parametricLengths( Cell( Dart( 30 ), 2 ) ), Eigen::Vector2d( 2.0, 0.3 ), 1e-10 ) );
}

TEST_CASE( "Multipatch parametric atlas delegates parametric starts to constituent patches" )
{
    const auto source_topo_0 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto line_topo_0 = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto tp_topo_0 = std::make_shared<const TPCombinatorialMap>( source_topo_0, line_topo_0 );
    const auto source_0 = std::make_shared<const ParametricAtlas1d>( source_topo_0, Eigen::Vector2d( 1.0, 2.0 ) );
    const auto line_0 = std::make_shared<const ParametricAtlas1d>( line_topo_0, Eigen::Vector3d( 0.5, 0.25, 0.125 ) );
    const auto tp_0 = std::make_shared<const TPParametricAtlas>( tp_topo_0, source_0, line_0 );

    const auto source_topo_1 = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto line_topo_1 = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto tp_topo_1 = std::make_shared<const TPCombinatorialMap>( source_topo_1, line_topo_1 );
    const auto source_1 = std::make_shared<const ParametricAtlas1d>( source_topo_1, Eigen::Vector2d( 4.0, 8.0 ) );
    const auto line_1 = std::make_shared<const ParametricAtlas1d>( line_topo_1, Eigen::Vector3d( 0.75, 1.5, 3.0 ) );
    const auto tp_1 = std::make_shared<const TPParametricAtlas>( tp_topo_1, source_1, line_1 );

    const auto mp_cmap = std::make_shared<const MultiPatchCombinatorialMap>(
        std::vector<std::shared_ptr<const TPCombinatorialMap>>{ tp_topo_0, tp_topo_1 },
        MultiPatchCombinatorialMap::InternalConnectionsMap{} );
    const MultiPatchParametricAtlas mp_atlas( mp_cmap, { tp_0, tp_1 } );

    const Cell local_patch_0_cell = topCellAt( *tp_topo_0, { 1, 2 } );
    const Cell global_patch_0_cell( mp_cmap->toGlobalDart( 0, local_patch_0_cell.dart() ), 2 );
    CHECK( util::equals( mp_atlas.parametricStarts( global_patch_0_cell ),
                         tp_0->parametricStarts( local_patch_0_cell ),
                         1e-12 ) );
    CHECK( util::equals( mp_atlas.parametricStarts( global_patch_0_cell ), Eigen::Vector2d( 1.0, 0.75 ), 1e-12 ) );

    const Cell local_patch_1_cell = topCellAt( *tp_topo_1, { 1, 1 } );
    const Cell global_patch_1_cell( mp_cmap->toGlobalDart( 1, local_patch_1_cell.dart() ), 2 );
    CHECK( util::equals( mp_atlas.parametricStarts( global_patch_1_cell ),
                         tp_1->parametricStarts( local_patch_1_cell ),
                         1e-12 ) );
    CHECK( util::equals( mp_atlas.parametricStarts( global_patch_1_cell ), Eigen::Vector2d( 4.0, 0.75 ), 1e-12 ) );
}

TEST_CASE( "corners" )
{
    const auto source_topo = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto line_topo = std::make_shared<const CombinatorialMap1d>( 2 );
    const auto tp_topo_2d = std::make_shared<const TPCombinatorialMap>( source_topo, line_topo );
    const auto tp_topo_3d = std::make_shared<const TPCombinatorialMap>( tp_topo_2d, line_topo );

    const auto source = std::make_shared<const ParametricAtlas1d>( source_topo, Eigen::Vector3d( 1, 2, 3 ) );
    const auto line = std::make_shared<const ParametricAtlas1d>( line_topo, Eigen::Vector2d( 0.1, 0.2 ) );
    const auto tp_2d = std::make_shared<const TPParametricAtlas>( tp_topo_2d, source, line );
    const TPParametricAtlas tp_3d( tp_topo_3d, tp_2d, line );


    const auto test_corners = []( const TPParametricAtlas& tp ) {
        const size_t dim = tp.cmap().dim();
        for( size_t i = 0; i < dim; i++ )
        {
            const auto corner_cells = cornerCells( tp, i );
            SmallVector< BaryCoordIsZeroVec, 12 > vecs;
            for( const auto& item : corner_cells )
            {
                const auto temp = parentDomainBoundary( tp, item );
                CHECK( std::find( vecs.begin(), vecs.end(), temp ) );
                vecs.push_back( temp );
                CHECK( tp.cmap().dim() - i == std::transform_reduce(
                                                  vecs.back().begin(),
                                                  vecs.back().end(),
                                                  size_t{ 0 },
                                                  []( const size_t accum, const size_t b ) { return accum + b; },
                                                  []( const bool b ) { return b ? 1 : 0; } ) );
            }
            CHECK( vecs.size() == ( i == 0 ? 1 : dim ) * pow( 2, dim - i ) );
        }
    };

    test_corners( *tp_2d );
    test_corners( tp_3d );
}
