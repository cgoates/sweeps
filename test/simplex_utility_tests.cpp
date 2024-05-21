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
