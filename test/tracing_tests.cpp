#include <catch2/catch_test_macros.hpp>
#include <Tracing.hpp>
#include <Logging.hpp>
#include <SimplexUtilities.hpp>
#include <CombinatorialMap.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>
#include <Laplace.hpp>
#include <MeshInput.hpp>
#include <CommonUtils.hpp>

namespace topology
{
class SingleTriangleCMap : public CombinatorialMap
{
    public:
    SingleTriangleCMap() {}
    virtual std::optional<Dart> phi( const int i, const Dart& d ) const override
    {
        if( std::abs( i ) != 1 ) return std::nullopt;
        return Dart( ( d.id() + 3 + i ) % 3 );
    }
    virtual Dart::IndexType maxDartId() const override { return 2; };
    virtual uint dim() const override { return 2; }
    virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override
    {
        for( const Dart::IndexType did : { 0, 1, 2 } ) if( not callback( Dart( did ) ) ) return false;
        return true;
    }
    virtual bool iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const override
    {
        switch( cell_dim )
        {
            case 2:
            return callback( Face( Dart( 0 ) ) );
            default:
            return iterateDartsWhile( [&]( const Dart& d ){ return callback( topology::Cell( d, cell_dim ) ); } );
        }
    }

    virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override
    {
        if( cell_dim != 0 ) return std::nullopt;
        else
        {
            return []( const Vertex& v ){
                return v.dart().id();
            };
        }
    }
};
}

TEST_CASE( "Test line ray intersection 0", "" )
{
    const Eigen::Vector2d line_point_0( 0.0, 1.0 );
    const Eigen::Vector2d line_point_1( 1.0, 0.0 );
    const Eigen::Vector2d ray_point( 0.0, 0.5 );
    const Eigen::Vector2d ray_dir( 0.1, 0.1 );

    const Segment<2> line{ line_point_0, line_point_1 };
    const Ray<2> ray{ ray_point, ray_dir };

    const auto intersection = reparam::intersectionOf( ray, line );

    REQUIRE( intersection.has_value() );
    REQUIRE( util::equals( intersection.value(), Eigen::Vector2d( 0.25, 0.75 ), 1e-8 ) );
}

TEST_CASE( "Test line ray intersection 1", "" )
{
    const Eigen::Vector2d line_point_1( 0.0, 1.0 );
    const Eigen::Vector2d line_point_0( 1.0, 0.0 );
    const Eigen::Vector2d ray_point( 0.0, 0.5 );
    const Eigen::Vector2d ray_dir( 0.1, 0.1 );

    const Segment<2> line{ line_point_0, line_point_1 };
    const Ray<2> ray{ ray_point, ray_dir };

    const auto intersection = reparam::intersectionOf( ray, line );

    REQUIRE( intersection.has_value() );
    REQUIRE( util::equals( intersection.value(), Eigen::Vector2d( 0.25, 0.75 ), 1e-8 ) );
}

TEST_CASE( "Test line ray no intersection", "" )
{
    const Eigen::Vector2d line_point_1( 0.0, 1.0 );
    const Eigen::Vector2d line_point_0( 1.0, 0.0 );
    const Eigen::Vector2d ray_point( -0.6, 0.5 );
    const Eigen::Vector2d ray_dir( 0.1, 0.1 );

    const Segment<2> line{ line_point_0, line_point_1 };
    const Ray<2> ray{ ray_point, ray_dir };

    const auto intersection = reparam::intersectionOf( ray, line );

    REQUIRE( not intersection.has_value() );
}

TEST_CASE( "Test line ray no intersection 2", "" )
{
    const Eigen::Vector2d line_point_1( 0.0, 1.0 );
    const Eigen::Vector2d line_point_0( 1.0, 0.0 );
    const Eigen::Vector2d ray_point( 0.0, 0.5 );
    const Eigen::Vector2d ray_dir( 0.1, -0.07 );

    const Segment<2> line{ line_point_0, line_point_1 };
    const Ray<2> ray{ ray_point, ray_dir };

    const auto intersection = reparam::intersectionOf( ray, line );

    REQUIRE( not intersection.has_value() );
}

TEST_CASE( "Test gradient tracing in triangle", "" )
{
    const std::array<Eigen::Vector3d, 3> tri(
        { Eigen::Vector3d( 1.0, 0.0, 0.0 ),
          Eigen::Vector3d( 0.0, 1.0, 0.0 ),
          Eigen::Vector3d( 0.0, 0.0, 1.0 ) } );
    const topology::SingleTriangleCMap map;
    const auto vertex_ids = indexingOrError( map, 0 );
    const auto positions = [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
        return tri.at( vertex_ids( v ) );
    };

    const topology::Edge e( topology::Dart( 0 ) );

    SECTION( "No intersection" )
    {
        const Eigen::Vector3d field_values( 0.0, 1.1, 0.5 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<topology::Edge, double>> intersection =
            reparam::traceGradientOnTri( map, positions, e, edge_barycentric_coord, field_values );
        REQUIRE( not intersection.has_value() );
    }

    SECTION( "Forward intersection" )
    {
        const Eigen::Vector3d field_values( 0.3, 1.1, 1.9 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<topology::Edge, double>> intersection =
            reparam::traceGradientOnTri( map, positions, e, edge_barycentric_coord, field_values );
        REQUIRE( intersection.has_value() );
        REQUIRE( intersection.value().first.dart().id() == 1 );
        REQUIRE( util::equals( intersection.value().second, 1.0 - edge_barycentric_coord, 1e-8 ) );
    }

    SECTION( "Backward intersection" )
    {
        const Eigen::Vector3d field_values( 1.1, 0.3, 1.9 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<topology::Edge, double>> intersection =
            reparam::traceGradientOnTri( map, positions, e, edge_barycentric_coord, field_values );
        REQUIRE( intersection.has_value() );
        REQUIRE( intersection.value().first.dart().id() == 2 );
        REQUIRE( util::equals( intersection.value().second, 1.0 - edge_barycentric_coord, 1e-8 ) );
    }

    // std::cout << "Intersection: {" << intersection.value().first << ", " << intersection.value().second << std::endl;
}