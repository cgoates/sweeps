#include <catch2/catch_test_macros.hpp>
#include <Tracing.hpp>
#include <Logging.hpp>

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
