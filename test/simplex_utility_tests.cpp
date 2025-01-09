#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <SimplexUtilities.hpp>
#include <Logging.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CommonUtils.hpp>

TEST_CASE( "Dihedral angle cotangent", "" )
{
    SimplicialComplex simplicial_complex = { { { 0, 1, 2, 3 } }, { { 0, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 } } };

    const auto test_angle =
        [&]( const double angle,
             const std::function<void( const topology::TetMeshCombinatorialMap&, const topology::Edge&, const std::vector<Normal>& )>&
                 test ) {
            simplicial_complex.points.push_back( { std::cos( angle ), std::sin( angle ), 0 } );

            topology::TetMeshCombinatorialMap map( simplicial_complex );

            const auto normals = faceNormals( map );

            const auto vertex_ids = indexingOrError( map, 0 );

            iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
                const auto vid1 = vertex_ids( topology::Vertex( e.dart() ) );
                const auto vid2 = vertex_ids( topology::Vertex( phi( map, 1, e.dart() ).value() ) );
                if( ( vid1 == 0 and vid2 == 1 ) or ( vid1 == 1 and vid2 == 0 ) )
                {
                    test( map, e, normals );
                    return false;
                }
                return true;
            } );
        };

    SECTION( "90 degrees" )
    {
        test_angle( std::numbers::pi / 2,
                    [&]( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals ) {
                        REQUIRE( util::equals( dihedralCotangent( map, e, normals ), 0, 1e-5 ) );
                    } );
    }

    SECTION( "45 degrees" )
    {
        test_angle( std::numbers::pi / 4,
                    [&]( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals ) {
                        REQUIRE( util::equals( dihedralCotangent( map, e, normals ), 1.0, 1e-5 ) );
                    } );
    }

    SECTION( "30 degrees" )
    {
        test_angle(
            std::numbers::pi / 6,
            [&]( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals ) {
                REQUIRE( util::equals( dihedralCotangent( map, e, normals ), 1.0 / std::tan( std::numbers::pi / 6 ), 1e-5 ) );
            } );
    }

    SECTION( "120 degrees" )
    {
        test_angle( 2 * std::numbers::pi / 3,
                    [&]( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals ) {
                        REQUIRE( util::equals(
                            dihedralCotangent( map, e, normals ), 1.0 / std::tan( 2 * std::numbers::pi / 3 ), 1e-5 ) );
                    } );
    }
}

TEST_CASE( "2d triangle inverse with canonical triangle" )
{
    const auto has_specific_value = []( const std::optional<Eigen::Vector3d>& inverse, const Eigen::Vector3d& expected ){
        CHECK( inverse.has_value() );
        if( inverse.has_value() ) CHECK( util::equals( inverse.value(), expected, 1e-12 ) );
    };
    SECTION( "Canonical triangle" )
    {
        const Triangle<2> tri{ { 0, 0 }, {1, 0}, {0, 1} };
        has_specific_value( invertTriangleMap( tri, {0.5, 0.5} ), {0.0, 0.5, 0.5} );
        has_specific_value( invertTriangleMap( tri, {0.0, 0.0} ), {1.0, 0.0, 0.0} );
        has_specific_value( invertTriangleMap( tri, {0.33, 0.33} ), {0.34, 0.33, 0.33} );
        has_specific_value( invertTriangleMap( tri, {0.21, 0.69} ), {0.1, 0.21, 0.69} );
        CHECK( not invertTriangleMap( tri, {1,1} ).has_value() );
        CHECK( not invertTriangleMap( tri, {-1, 0.2} ).has_value() );
    }

    SECTION( "" )
    {
        const Triangle<2> tri{ { 0.23, 2 }, {1.9, 0.7}, {2.3, 1.8} };
        has_specific_value( invertTriangleMap( tri, {2.1, 1.25} ), {0.0, 0.5, 0.5} );
        has_specific_value( invertTriangleMap( tri, {0.23, 2.0} ), {1.0, 0.0, 0.0} );
        has_specific_value( invertTriangleMap( tri, {1.4642, 1.505} ), {0.34, 0.33, 0.33} );
        has_specific_value( invertTriangleMap( tri, {2.009, 1.589} ), {0.1, 0.21, 0.69} );

        CHECK( not invertTriangleMap( tri, {1,1} ).has_value() );
        CHECK( not invertTriangleMap( tri, {-1, 0.2} ).has_value() );
        CHECK( not invertTriangleMap( tri, {2.2, 1.25} ).has_value() );
    }
}

TEST_CASE( "3d triangle inverse" )
{
    const auto has_specific_value = []( const std::optional<Eigen::Vector3d>& inverse,
                                        const Eigen::Vector3d& expected ) {
        CHECK( inverse.has_value() );
        if( inverse.has_value() ) CHECK( util::equals( inverse.value(), expected, 1e-12 ) );
    };
    SECTION( "Canonical triangle" )
    {
        const Triangle<3> tri{ { 0, 0, 1 }, { 1, 0, 0 }, { 0, 1, 0 } };
        has_specific_value( invertTriangleMap( tri, { 0.5, 0.5, 0.0 } ), { 0.0, 0.5, 0.5 } );
        has_specific_value( invertTriangleMap( tri, { 0.0, 0.0, 1.0 } ), { 1.0, 0.0, 0.0 } );
        has_specific_value( invertTriangleMap( tri, { 0.3, 0.4, 0.3 } ), { 0.3, 0.3, 0.4 } );
    }
    SECTION( "Random triangle" )
    {
        const Triangle<3> tri{ { 0.23, 2, 0.7 }, { 1.9, 0.7, 6.2 }, { 2.3, 1.8, -10.0 } };
        has_specific_value( invertTriangleMap( tri, { 1.065, 1.35, 3.45 } ), { 0.5, 0.5, 0.0 } );
        has_specific_value( invertTriangleMap( tri, { 0.23, 2.0, 0.7 } ), { 1.0, 0.0, 0.0 } );
        has_specific_value( invertTriangleMap( tri, { 1.766, 1.51, -3.0 } ), { 0.2, 0.3, 0.5 } );
    }
}

TEST_CASE( "Triangle closest point" )
{
    const auto in_triangle_with_value = []( const std::pair<Eigen::Vector3d, std::optional<double>>& closest, const Eigen::Vector3d& expected ){
        CHECK( not closest.second.has_value() );
        CHECK( util::equals( closest.first, expected, 1e-12 ) );
    };
    const auto out_of_triangle_with_value = []( const std::pair<Eigen::Vector3d, std::optional<double>>& closest, const Eigen::Vector3d& expected ){
        CHECK( closest.second.has_value() );
        CHECK( util::equals( closest.first, expected, 1e-12 ) );
    };
    SECTION( "Canonical triangle" )
    {
        const Triangle<3> tri{ { 0, 0, 1 }, { 1, 0, 0 }, { 0, 1, 0 } };
        in_triangle_with_value( invertTriangleMapOrClosestPoint( tri, { 0.3, 0.4, 0.3 } ), { 0.3, 0.3, 0.4 } );
        out_of_triangle_with_value( invertTriangleMapOrClosestPoint( tri, { 0, 0, 1.1 } ), { 1, 0, 0 } );
        out_of_triangle_with_value( invertTriangleMapOrClosestPoint( tri, { 0.5, 0.5, 0.5 } ), { 1.0/3.0, 1.0/3.0, 1.0/3.0 } );
        out_of_triangle_with_value( invertTriangleMapOrClosestPoint( tri, { 0.5, 0.5, -0.1 } ), { 0.0, 0.5, 0.5 } );
    }
}