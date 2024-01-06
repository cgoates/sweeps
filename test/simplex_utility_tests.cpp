#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <SimplexUtilities.hpp>

TEST_CASE( "Dihedral angle cotangent", "[single-file]" ) {
    SweepInput sweep_input = { {{0, 1, 2, 3}}, { {0, 0, 0}, {0, 0, 1}, {1, 0, 0} }, {}, {} };

    const auto test_angle = [&]( const double angle,
                                 const std::function<void( const cgogn::CMap3&,
                                                           const cgogn::CMap3::Edge&,
                                                           const std::vector<Normal>& )>& test ) {
        sweep_input.points.push_back( {std::cos( angle ), std::sin( angle ), 0} );

        cgogn::CMap3 map;
        SimplexUtilities::mapFromInput( sweep_input, map );

        const auto normals = SimplexUtilities::faceNormals( map );

        cgogn::foreach_cell( map, [&]( cgogn::CMap3::Edge e ) {
            const auto vid1 = cgogn::index_of( map, cgogn::CMap3::Vertex( e.dart_ ) );
            const auto vid2 = cgogn::index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, e.dart_ ) ) );
            if( ( vid1 == 0 and vid2 == 1 ) or ( vid1 == 1 and vid2 == 0 ) )
            {
                test( map, e, normals );
                return false;
            }
            return true;
        } );
    };

    SECTION( "90 degrees" ){
        test_angle( std::numbers::pi/2, [&]( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e, const std::vector<Normal>& normals ) {
            REQUIRE( SimplexUtilities::dihedralCotangent( map, e, normals ) < 1e-5 );
        } );
    }

    SECTION( "45 degrees" ){
        test_angle( std::numbers::pi/4, [&]( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e, const std::vector<Normal>& normals ) {
            REQUIRE( std::abs( SimplexUtilities::dihedralCotangent( map, e, normals ) - 1.0 ) < 1e-5 );
        } );
    }

    SECTION( "30 degrees" ){
        test_angle( std::numbers::pi/6, [&]( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e, const std::vector<Normal>& normals ) {
            REQUIRE( std::abs( SimplexUtilities::dihedralCotangent( map, e, normals ) - 1.0/std::tan( std::numbers::pi/6 ) ) < 1e-5 );
        } );
    }

    SECTION( "120 degrees" ){
        test_angle( 2*std::numbers::pi/3, [&]( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e, const std::vector<Normal>& normals ) {
            REQUIRE( std::abs( SimplexUtilities::dihedralCotangent( map, e, normals ) - 1.0/std::tan( 2*std::numbers::pi/3 ) ) < 1e-5 );
        } );
    }
}

