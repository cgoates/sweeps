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
#include <AbaqusInput.hpp>

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

    virtual VertexId vertexId( const Vertex& v ) const override
    {
        return v.dart().id();
    }
};
}

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
    const std::array<Eigen::Vector3d, 3> tri(
        { Eigen::Vector3d( 1.0, 0.0, 0.0 ),
          Eigen::Vector3d( 0.0, 1.0, 0.0 ),
          Eigen::Vector3d( 0.0, 0.0, 1.0 ) } );
    const topology::SingleTriangleCMap map;
    const auto positions = [&]( const topology::Vertex& v ) -> const Eigen::Vector3d& {
        return tri.at( map.vertexId( v ).id() );
    };

    const topology::Edge e( topology::Dart( 0 ) );

    SECTION( "No intersection" )
    {
        const Eigen::Vector3d field_values( 0.0, 1.1, 0.5 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<topology::Edge, double>> intersection =
            traceGradientOnTri( map, positions, e, edge_barycentric_coord, field_values );
        REQUIRE( not intersection.has_value() );
    }

    SECTION( "Forward intersection" )
    {
        const Eigen::Vector3d field_values( 0.3, 1.1, 1.9 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<topology::Edge, double>> intersection =
            traceGradientOnTri( map, positions, e, edge_barycentric_coord, field_values );
        REQUIRE( intersection.has_value() );
        REQUIRE( intersection.value().first.dart().id() == 1 );
        REQUIRE( equals( intersection.value().second, 1.0 - edge_barycentric_coord, 1e-8 ) );
    }

    SECTION( "Backward intersection" )
    {
        const Eigen::Vector3d field_values( 1.1, 0.3, 1.9 );
        const double edge_barycentric_coord = 0.7;
        const std::optional<std::pair<topology::Edge, double>> intersection =
            traceGradientOnTri( map, positions, e, edge_barycentric_coord, field_values );
        REQUIRE( intersection.has_value() );
        REQUIRE( intersection.value().first.dart().id() == 2 );
        REQUIRE( equals( intersection.value().second, 1.0 - edge_barycentric_coord, 1e-8 ) );
    }

    // std::cout << "Intersection: {" << intersection.value().first << ", " << intersection.value().second << std::endl;
}

TEST_CASE( "Tracing from all the cells in the macaroni", "[slow]")
{
    const SweepInput sweep_input =
        io::loadINPFile( SRC_HOME "/test/data/macaroni.inp", "Surface3", "Surface4" );

    const topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    std::cout << "Solving laplace\n";
    const Eigen::VectorXd ans =
        solveLaplaceSparse( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

    const topology::CombinatorialMapBoundary bdry( map );

    const auto keep_face_sides = [&]( const topology::Face& f ) {
        return ( not sweep_input.zero_bcs.at( bdry.vertexId( topology::Vertex( f.dart() ) ).id() ) or
                    not sweep_input.zero_bcs.at(
                        bdry.vertexId( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ).id() ) or
                    not sweep_input.zero_bcs.at(
                        bdry.vertexId( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ).id() ) ) and
                ( not sweep_input.one_bcs.at( bdry.vertexId( topology::Vertex( f.dart() ) ).id() ) or
                    not sweep_input.one_bcs.at(
                        bdry.vertexId( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ).id() ) or
                    not sweep_input.one_bcs.at(
                        bdry.vertexId( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ).id() ) );
    };

    const auto keep_face_base = [&]( const topology::Face& f ) {
        return sweep_input.zero_bcs.at( bdry.vertexId( topology::Vertex( f.dart() ) ).id() ) and
                sweep_input.zero_bcs.at( bdry.vertexId( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ).id() ) and
                sweep_input.zero_bcs.at( bdry.vertexId( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ).id() );
    };

    const auto keep_face_target = [&]( const topology::Face& f ) {
        return sweep_input.one_bcs.at( bdry.vertexId( topology::Vertex( f.dart() ) ).id() ) and
                sweep_input.one_bcs.at(
                    bdry.vertexId( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ).id() ) and
                sweep_input.one_bcs.at(
                    bdry.vertexId( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ).id() );
    };

    const topology::CombinatorialMapRestriction sides( bdry, keep_face_sides );
    const topology::CombinatorialMapRestriction base( bdry, keep_face_base );
    const topology::CombinatorialMapRestriction target( bdry, keep_face_target );

    const Eigen::Matrix3Xd grad = gradientsWithBoundaryCorrection( map, sides, ans, normals );

    const auto vertex_positions = [&sweep_input]( const topology::CombinatorialMap& map ){
        return [&sweep_input, &map]( const topology::Vertex& v ) -> const Eigen::Vector3d& {
            return sweep_input.mesh.points.at( map.vertexId( v ).id() );
        };
    };

    iterateCellsWhile( map, 2, [&]( const topology::Face& start_face ) {
        if( onBoundary( map, start_face.dart() ) ) return true;
        CHECK_NOTHROW( traceField( map, start_face, centroid( map, start_face ), grad, normals ) );
        return true;
    } );
    iterateCellsWhile( base, 2, [&]( const topology::Face& start_face ) {
        CHECK_NOTHROW( traceField( map, start_face, centroid( map, start_face ), grad, normals ) );
        return true;
    } );

    iterateCellsWhile( sides, 1, [&]( const topology::Edge& start_edge ) {
        const auto maybe_phi2 = phi( bdry, 2, start_edge.dart() );
        if( maybe_phi2.has_value() and keep_face_target( maybe_phi2.value() ) ) return true;
        CHECK_NOTHROW( traceBoundaryField( sides, start_edge, 0.5, ans, vertex_positions( sides ), false ) );
        return true;
    } );

    const auto bdry_positions = vertex_positions( bdry );
    topology::GlobalCellMarker traced_vertices( map, 0 );

    iterateCellsWhile( base, 0, [&]( const topology::Vertex& v ) {
        traced_vertices.mark( map, bdry.toUnderlyingCell( v ) );
        return true;
    } );

    // Tracing reverse as a parameterization test
    const auto reverse_ans = -1 * ans;
    const auto reverse_grad = -1 * grad;

    iterateCellsWhile( sides, 0, [&]( const topology::Vertex& v ) {
        if( traced_vertices.isMarked( bdry.toUnderlyingCell( v ) ) ) return true;
        traced_vertices.mark( map, bdry.toUnderlyingCell( v ) );
        CHECK_NOTHROW( traceBoundaryField( sides, v, 1.0, reverse_ans, bdry_positions, false ) );
        return true;
    } );

    const auto pos = vertex_positions( map );
    iterateCellsWhile( map, 0, [&]( const topology::Vertex& v ) {
        if( traced_vertices.isMarked( v ) ) return true;
        CHECK_NOTHROW( traceField( map, v, pos( v ), reverse_grad, normals, false ) );
        return true;
    } );
}
