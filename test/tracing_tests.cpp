#include <catch2/catch_test_macros.hpp>
#include <Tracing.hpp>
#include <Logging.hpp>
#include <SimplexUtilities.hpp>

TEST_CASE( "Test line ray intersection 0", "[single-file]" )
{
    const Eigen::Vector2d line_point_0( 0.0, 1.0 );
    const Eigen::Vector2d line_point_1( 1.0, 0.0 );
    const Eigen::Vector2d ray_point( 0.0, 0.5 );
    const Eigen::Vector2d ray_dir( 0.1, 0.1 );

    const Segment<2> line{ line_point_0, line_point_1 };
    const Ray<2> ray{ ray_point, ray_dir };

    const auto intersection = intersectionOf( ray, line );

    REQUIRE( intersection.has_value() );
    REQUIRE( equals( intersection.value(), Eigen::Vector2d( 0.25, 0.75 ), 1e-8 ) );
}

TEST_CASE( "Test line ray intersection 1", "[single-file]" )
{
    const Eigen::Vector2d line_point_1( 0.0, 1.0 );
    const Eigen::Vector2d line_point_0( 1.0, 0.0 );
    const Eigen::Vector2d ray_point( 0.0, 0.5 );
    const Eigen::Vector2d ray_dir( 0.1, 0.1 );

    const Segment<2> line{ line_point_0, line_point_1 };
    const Ray<2> ray{ ray_point, ray_dir };

    const auto intersection = intersectionOf( ray, line );

    REQUIRE( intersection.has_value() );
    REQUIRE( equals( intersection.value(), Eigen::Vector2d( 0.25, 0.75 ), 1e-8 ) );
}

TEST_CASE( "Test line ray no intersection", "[single-file]" )
{
    const Eigen::Vector2d line_point_1( 0.0, 1.0 );
    const Eigen::Vector2d line_point_0( 1.0, 0.0 );
    const Eigen::Vector2d ray_point( -0.6, 0.5 );
    const Eigen::Vector2d ray_dir( 0.1, 0.1 );

    const Segment<2> line{ line_point_0, line_point_1 };
    const Ray<2> ray{ ray_point, ray_dir };

    const auto intersection = intersectionOf( ray, line );

    REQUIRE( not intersection.has_value() );
}

TEST_CASE( "Test line ray no intersection 2", "[single-file]" )
{
    const Eigen::Vector2d line_point_1( 0.0, 1.0 );
    const Eigen::Vector2d line_point_0( 1.0, 0.0 );
    const Eigen::Vector2d ray_point( 0.0, 0.5 );
    const Eigen::Vector2d ray_dir( 0.1, -0.07 );

    const Segment<2> line{ line_point_0, line_point_1 };
    const Ray<2> ray{ ray_point, ray_dir };

    const auto intersection = intersectionOf( ray, line );

    REQUIRE( not intersection.has_value() );
}

TEST_CASE( "Test gradient tracing in triangle", "[single-file]" )
{
    const Eigen::Vector3d v1( 1.0, 0.0, 0.0 );
    const Eigen::Vector3d v2( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d v3( 0.0, 0.0, 1.0 );
    const Triangle<3> tri( { v1, v2, v3 } );

    SECTION( "No intersection" )
    {
        const Eigen::Vector3d field_values( 0.0, 1.1, 0.5 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<bool, double>> intersection =
            traceGradientOnTri( tri, edge_barycentric_coord, field_values );
        REQUIRE( not intersection.has_value() );
    }

    SECTION( "Forward intersection" )
    {
        const Eigen::Vector3d field_values( 0.3, 1.1, 1.9 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<bool, double>> intersection =
            traceGradientOnTri( tri, edge_barycentric_coord, field_values );
        REQUIRE( intersection.has_value() );
        REQUIRE( intersection.value().first );
        REQUIRE( equals( intersection.value().second, 1.0 - edge_barycentric_coord, 1e-8 ) );
    }

    SECTION( "Backward intersection" )
    {
        const Eigen::Vector3d field_values( 1.1, 0.3, 1.9 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<bool, double>> intersection =
            traceGradientOnTri( tri, edge_barycentric_coord, field_values );
        REQUIRE( intersection.has_value() );
        REQUIRE( not intersection.value().first );
        REQUIRE( equals( intersection.value().second, 1.0 - edge_barycentric_coord, 1e-8 ) );
    }

    // std::cout << "Intersection: {" << intersection.value().first << ", " << intersection.value().second << std::endl;
}
