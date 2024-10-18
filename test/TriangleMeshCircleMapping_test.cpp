#include <catch2/catch_test_macros.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <TriMeshCombinatorialMap.hpp>
#include <TriangleParametricAtlas.hpp>
#include <TriangleMeshCircleMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParentDomain.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

TEST_CASE( "Simple circle evaluations" )
{
    // This but fit boundary edges to a circle:
    //   *
    //  /|\.
    // *-*-*
    //  \|/
    //   *
    const SimplicialComplex mesh{
        { { 0, 1, 4 }, { 0, 4, 2 }, { 4, 1, 3 }, { 4, 3, 2 } },
        { { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.2 }, { 0.0, 0.0, 0.0 } } };

    const auto cmap = std::make_shared<const topology::TriMeshCombinatorialMap>( mesh );
    const auto vert_id = indexingOrError( *cmap, 0 );
    const auto param = std::make_shared<const param::TriangleParametricAtlas>( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) ).head<2>();
    };
    const mapping::TriangleMeshCircleMapping geom_mapping( param, vert_positions );

    const param::ParentDomain pd = param::simplexDomain( 2 );
    iterateCellsWhile( *cmap, 2, [&]( const topology::Face& f ) {
        iterateAdjacentCells( *cmap, f, 1, [&]( const topology::Edge& e ) {
            const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( *param, e ) );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            if( boundaryAdjacent( *cmap, e ) )
                CHECK( util::equals( pos.norm(), 1.0, 1e-12 ) );
            else
                CHECK( util::equals( pos.norm(), 0.5, 1e-12 ) );
            return true;
        } );

        iterateAdjacentCells( *cmap, f, 0, [&]( const topology::Vertex& v ) {
            const param::ParentPoint ppt = param->parentPoint( v );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const Eigen::Vector2d expected = vert_positions( v );
            CHECK( util::equals( pos, expected, 1e-12 ) );
            return true;
        } );

        const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( *param, f ) );
        const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
        CHECK( util::equals( pos.norm(), 2.0/3.0, 1e-12 ) );
        return true;
    } );
}

TEST_CASE( "Circle evaluation with non-boundary-adjacent cells" )
{
    // An equilateral triangle within a circle
    std::vector<Eigen::Vector3d> points;
    for( const double angle : { 0.0, 2.0 * std::numbers::pi / 3.0, 4.0 * std::numbers::pi / 3.0 } )
    {
        points.emplace_back( 0.5 * std::cos( angle ), 0.5 * std::sin( angle ), 0.0 );
        points.emplace_back( std::cos( angle ), std::sin( angle ), 0.0 );
    }
    const SimplicialComplex mesh{
        { {0, 2, 4}, {0, 1, 2}, {2, 1, 3}, {2, 3, 4}, {4, 3, 5}, {4, 5, 0}, {0, 5, 1} },
        points
    };

        const auto cmap = std::make_shared<const topology::TriMeshCombinatorialMap>( mesh );
    const auto vert_id = indexingOrError( *cmap, 0 );
    const auto param = std::make_shared<const param::TriangleParametricAtlas>( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) ).head<2>();
    };
    const mapping::TriangleMeshCircleMapping geom_mapping( param, vert_positions );

    const param::ParentDomain pd = param::simplexDomain( 2 );
    iterateCellsWhile( *cmap, 2, [&]( const topology::Face& f ) {
        iterateAdjacentCells( *cmap, f, 1, [&]( const topology::Edge& e ) {
            const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( *param, e ) );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            if( boundaryAdjacent( *cmap, e ) )
                CHECK( util::equals( pos.norm(), 1.0, 1e-12 ) );
            return true;
        } );

        iterateAdjacentCells( *cmap, f, 0, [&]( const topology::Vertex& v ) {
            const param::ParentPoint ppt = param->parentPoint( v );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const Eigen::Vector2d expected = vert_positions( v );
            CHECK( util::equals( pos, expected, 1e-12 ) );
            return true;
        } );
        return true;
    } );

    const param::ParentPoint ppt = pointOnBoundary( pd, {false, false, false} );
    const Eigen::Vector2d pos = geom_mapping.evaluate( topology::Face( topology::Dart( 0 ) ), ppt );
    CHECK( util::equals( pos, Eigen::Vector2d( 0.0, 0.0 ), 1e-12 ) );
}

TEST_CASE( "Circle mapping and inverse round trip test" )
{
    // An equilateral triangle within a circle
    std::vector<Eigen::Vector3d> points;
    for( const double angle : { 0.0, 2.0 * std::numbers::pi / 3.0, 4.0 * std::numbers::pi / 3.0 } )
    {
        points.emplace_back( 0.5 * std::cos( angle ), 0.5 * std::sin( angle ), 0.0 );
        points.emplace_back( std::cos( angle ), std::sin( angle ), 0.0 );
    }
    const SimplicialComplex mesh{
        { {0, 2, 4}, {0, 1, 2}, {2, 1, 3}, {2, 3, 4}, {4, 3, 5}, {4, 5, 0}, {0, 5, 1} },
        points
    };

    const auto cmap = std::make_shared<const topology::TriMeshCombinatorialMap>( mesh );
    const auto vert_id = indexingOrError( *cmap, 0 );
    const auto param = std::make_shared<const param::TriangleParametricAtlas>( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) ).head<2>();
    };
    const mapping::TriangleMeshCircleMapping geom_mapping( param, vert_positions );

    const param::ParentDomain pd = param::simplexDomain( 2 );
    iterateCellsWhile( *cmap, 2, [&]( const topology::Face& f ) {
        iterateAdjacentCells( *cmap, f, 1, [&]( const topology::Edge& e ) {
            const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( *param, e ) );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const auto ppt_rt = geom_mapping.maybeInverse( f, pos );
            CHECK( ppt_rt.has_value() );
            if( ppt_rt.has_value() )
            {
                CHECK( ppt_rt.value().mBaryCoordIsZero == ppt.mBaryCoordIsZero );
                CHECK( util::equals( ppt_rt->mPoint, ppt.mPoint, 1e-10 ) );
            }
            return true;
        } );

        iterateAdjacentCells( *cmap, f, 0, [&]( const topology::Vertex& v ) {
            const param::ParentPoint ppt = param->parentPoint( v );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const auto ppt_rt = geom_mapping.maybeInverse( f, pos );
            CHECK( ppt_rt.has_value() );
            if( ppt_rt.has_value() )
            {
                CHECK( ppt_rt.value().mBaryCoordIsZero == ppt.mBaryCoordIsZero );
                CHECK( util::equals( ppt_rt->mPoint, ppt.mPoint, 1e-10 ) );

            }
            return true;
        } );

        const param::ParentPoint ppt = pointOnBoundary( pd, {false, false, false} );
        const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
        const auto ppt_rt = geom_mapping.maybeInverse( f, pos );
        CHECK( ppt_rt.has_value() );
        if( ppt_rt.has_value() )
        {
            CHECK( ppt_rt.value().mBaryCoordIsZero == ppt.mBaryCoordIsZero );
            CHECK( util::equals( ppt_rt->mPoint, ppt.mPoint, 1e-10 ) );
        }
        return true;
    } );
}

TEST_CASE( "Mesh wide circle map inverse test" )
{
    // This but fit boundary edges to a circle:
    //   *
    //  /|\.
    // *-*-*
    //  \|/
    //   *
    const SimplicialComplex mesh{
        { { 0, 1, 4 }, { 0, 4, 2 }, { 4, 1, 3 }, { 4, 3, 2 } },
        { { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } } };

    const auto cmap = std::make_shared<const topology::TriMeshCombinatorialMap>( mesh );
    const auto vert_id = indexingOrError( *cmap, 0 );
    const auto param = std::make_shared<const param::TriangleParametricAtlas>( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) ).head<2>();
    };
    const mapping::TriangleMeshCircleMapping geom_mapping( param, vert_positions );

    const auto face_ids = indexingOrError( *cmap, 2 );

    const auto test_inverse = [&]( const Eigen::Vector2d& spatial_pt,
                                   const topology::Face& expected_f,
                                   const Eigen::Vector2d& expected_ppt ) {
        const auto inverse_pr = geom_mapping.maybeInverse( spatial_pt );
        CHECK( inverse_pr.has_value() );
        if( inverse_pr.has_value() )
        {
            CHECK( face_ids( expected_f ) == face_ids( inverse_pr->first ) );
            std::cout << expected_ppt.transpose() << " | " << inverse_pr->second.mPoint.transpose() << std::endl;
            CHECK( util::equals( expected_ppt, inverse_pr->second.mPoint, 1e-5 ) );
        }
    };
    const auto test_no_inverse = [&]( const Eigen::Vector2d& spatial_pt ){
        CHECK( not geom_mapping.maybeInverse( spatial_pt ).has_value() );
    };

    test_no_inverse( {4.4, 0} );
    test_inverse( {-0.25, -0.25}, topology::Face( topology::Dart( 0 ) ), {0.176776695, 0.6464466094} );
    test_inverse( { 0.25, -0.25}, topology::Face( topology::Dart( 6 ) ), {0.176776695, 0.176776695} );
    test_inverse( {-0.25,  0.25}, topology::Face( topology::Dart( 3 ) ), {0.6464466094, 0.176776695} );
    test_inverse( { 0.25,  0.25}, topology::Face( topology::Dart( 9 ) ), {0.176776695, 0.176776695} );
    test_no_inverse( {0.9, 0.5} );
}

TEST_CASE( "Simple circle evaluations two boundary edges" )
{
    // This but fit boundary edges to a circle:
    //   *
    //  / \.
    // *---*
    //  \ /
    //   *
    const SimplicialComplex mesh{
        { { 0, 1, 3 }, { 0, 3, 2 } },
        { { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 } } };

    const auto cmap = std::make_shared<const topology::TriMeshCombinatorialMap>( mesh );
    const auto vert_id = indexingOrError( *cmap, 0 );
    const auto param = std::make_shared<const param::TriangleParametricAtlas>( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) ).head<2>();
    };
    const mapping::TriangleMeshCircleMapping geom_mapping( param, vert_positions );

    const param::ParentDomain pd = param::simplexDomain( 2 );
    iterateCellsWhile( *cmap, 2, [&]( const topology::Face& f ) {
        iterateAdjacentCells( *cmap, f, 1, [&]( const topology::Edge& e ) {
            const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( *param, e ) );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            if( boundaryAdjacent( *cmap, e ) )
                CHECK( util::equals( pos.norm(), 1.0, 1e-12 ) );
            else
                CHECK( util::equals( pos.norm(), 0.0, 1e-12 ) );
            return true;
        } );

        iterateAdjacentCells( *cmap, f, 0, [&]( const topology::Vertex& v ) {
            const param::ParentPoint ppt = param->parentPoint( v );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const Eigen::Vector2d expected = vert_positions( v );
            CHECK( util::equals( pos, expected, 1e-12 ) );
            return true;
        } );

        const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( *param, f ) );
        const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
        CHECK( util::equals( pos.norm(), 1.0/3.0, 1e-12 ) );
        return true;
    } );
}

TEST_CASE( "Circle mapping and inverse round trip test with two boundary edge tris" )
{
    // This but fit boundary edges to a circle:
    //   *
    //  / \.
    // *---*
    //  \ /
    //   *
    const SimplicialComplex mesh{
        { { 0, 1, 3 }, { 0, 3, 2 } },
        { { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 } } };

    const auto cmap = std::make_shared<const topology::TriMeshCombinatorialMap>( mesh );
    const auto vert_id = indexingOrError( *cmap, 0 );
    const auto param = std::make_shared<const param::TriangleParametricAtlas>( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) ).head<2>();
    };
    const mapping::TriangleMeshCircleMapping geom_mapping( param, vert_positions );

    const param::ParentDomain pd = param::simplexDomain( 2 );
    iterateCellsWhile( *cmap, 2, [&]( const topology::Face& f ) {
        iterateAdjacentCells( *cmap, f, 1, [&]( const topology::Edge& e ) {
            const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( *param, e ) );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const auto ppt_rt = geom_mapping.maybeInverse( f, pos );
            CHECK( ppt_rt.has_value() );
            if( ppt_rt.has_value() )
            {
                CHECK( ppt_rt.value().mBaryCoordIsZero == ppt.mBaryCoordIsZero );
                CHECK( util::equals( ppt_rt->mPoint, ppt.mPoint, 1e-10 ) );
            }
            return true;
        } );

        iterateAdjacentCells( *cmap, f, 0, [&]( const topology::Vertex& v ) {
            const param::ParentPoint ppt = param->parentPoint( v );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const auto ppt_rt = geom_mapping.maybeInverse( f, pos );
            CHECK( ppt_rt.has_value() );
            if( ppt_rt.has_value() )
            {
                CHECK( ppt_rt.value().mBaryCoordIsZero == ppt.mBaryCoordIsZero );
                CHECK( util::equals( ppt_rt->mPoint, ppt.mPoint, 1e-10 ) );

            }
            return true;
        } );

        const param::ParentPoint ppt = pointOnBoundary( pd, {false, false, false} );
        const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
        const auto ppt_rt = geom_mapping.maybeInverse( f, pos );
        CHECK( ppt_rt.has_value() );
        if( ppt_rt.has_value() )
        {
            CHECK( ppt_rt.value().mBaryCoordIsZero == ppt.mBaryCoordIsZero );
            CHECK( util::equals( ppt_rt->mPoint, ppt.mPoint, 1e-10 ) );
        }
        return true;
    } );
}