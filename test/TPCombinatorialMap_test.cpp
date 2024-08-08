#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>

using namespace topology;
TEST_CASE( "3x4 TP combinatorial map" )
{
    const auto source = std::make_shared<const CombinatorialMap1d>( 3 );
    const auto line = std::make_shared<const CombinatorialMap1d>( 4 );
    const TPCombinatorialMap tp( source, line );

    const std::array<Dart, 48> phi_ones(
        { 1,  2,  3,  0,  5,  6,  7,  4,  9,  10, 11, 8,  13, 14, 15, 12, 17, 18, 19, 16, 21, 22, 23, 20,
          25, 26, 27, 24, 29, 30, 31, 28, 33, 34, 35, 32, 37, 38, 39, 36, 41, 42, 43, 40, 45, 46, 47, 44 } );
    const std::array<Dart, 48> phi_minus_ones(
        { 3,  0,  1,  2,  7,  4,  5,  6,  11, 8,  9,  10, 15, 12, 13, 14, 19, 16, 17, 18, 23, 20, 21, 22,
          27, 24, 25, 26, 31, 28, 29, 30, 35, 32, 33, 34, 39, 36, 37, 38, 43, 40, 41, 42, 47, 44, 45, 46 } );

    const std::array<std::optional<Dart>, 48> phi_twos(
        { std::nullopt, 7,  12, std::nullopt, std::nullopt, 11, 16, 1,  std::nullopt, std::nullopt,
          20, 5,  2,  19, 24, std::nullopt, 6,  23, 28, 13, 10, std::nullopt, 32, 17, 14, 31, 36,
          std::nullopt, 18, 35, 40, 25, 22, std::nullopt, 44, 29, 26, 43, std::nullopt, std::nullopt,
          30, 47, std::nullopt, 37, 34, std::nullopt, std::nullopt, 41 } );

    iterateDartsWhile( tp, [&]( const Dart& d ){
        const auto phi1 = phi( tp, 1, d );
        const auto phi_1 = phi( tp, -1, d );
        const auto phi2 = phi( tp, 2, d );
        CHECK( phi1.has_value() );
        CHECK( phi_1.has_value() );
        if( phi1.has_value() ) CHECK( phi1.value() == phi_ones.at( d.id() ) );
        if( phi_1.has_value() ) CHECK( phi_1.value() == phi_minus_ones.at( d.id() ) );
        CHECK( phi2 == phi_twos.at( d.id() ) );
        return true;
    } );
}

TEST_CASE( "2x2 TP combinatorial map" )
{
    const auto two_elems_1d = std::make_shared<const CombinatorialMap1d>( 2 );
    const TPCombinatorialMap tp( two_elems_1d, two_elems_1d );

    const std::array<Dart, 16> phi_ones( { 1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12 } );
    const std::array<Dart, 16> phi_minus_ones( { 3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10, 15, 12, 13, 14 } );

    const std::array<std::optional<Dart>, 16> phi_twos(
        { std::nullopt, 7,  8, std::nullopt, std::nullopt, std::nullopt, 12, 1,
          2, 15, std::nullopt, std::nullopt, 6, std::nullopt, std::nullopt, 9 } );

    iterateDartsWhile( tp, [&]( const Dart& d ){
        const auto phi1 = phi( tp, 1, d );
        const auto phi_1 = phi( tp, -1, d );
        const auto phi2 = phi( tp, 2, d );
        CHECK( phi1.has_value() );
        CHECK( phi_1.has_value() );
        if( phi1.has_value() ) CHECK( phi1.value() == phi_ones.at( d.id() ) );
        if( phi_1.has_value() ) CHECK( phi_1.value() == phi_minus_ones.at( d.id() ) );
        CHECK( phi2 == phi_twos.at( d.id() ) );
        return true;
    } );
}

TEST_CASE( "3d TP combinatorial map" )
{
    const auto one_elem_1d = std::make_shared<const CombinatorialMap1d>( 1 );
    const auto tp_2d = std::make_shared<const TPCombinatorialMap>( one_elem_1d, one_elem_1d );
    const TPCombinatorialMap tp( tp_2d, one_elem_1d );

    const std::array<Dart, 24> phi_ones(
        { 6, 2, 3, 4, 1, 23, 12, 8, 9, 10, 7, 5, 18, 14, 15, 16, 13, 11, 0, 20, 21, 22, 19, 17 } );
    const std::array<Dart, 24> phi_minus_ones(
        { 18, 4, 1, 2, 3, 11, 0, 10, 7, 8, 9, 17, 6, 16, 13, 14, 15, 23, 12, 22, 19, 20, 21, 5 } );
    const std::array<Dart, 48> phi_twos(
        { 1, 0, 22, 5, 8, 3, 7, 6, 4, 11, 14, 9, 13, 12, 10, 17, 20, 15, 19, 18, 16, 23, 2, 21 } );

    iterateDartsWhile( tp, [&]( const Dart& d ){
        const auto phi1 = phi( tp, 1, d );
        const auto phi_1 = phi( tp, -1, d );
        const auto phi2 = phi( tp, 2, d );
        CHECK( phi1.has_value() );
        CHECK( phi_1.has_value() );
        CHECK( phi2.has_value() );
        if( phi1.has_value() ) CHECK( phi1.value() == phi_ones.at( d.id() ) );
        if( phi_1.has_value() ) CHECK( phi_1.value() == phi_minus_ones.at( d.id() ) );
        if( phi2.has_value() ) CHECK( phi2.value() == phi_twos.at( d.id() ) );
        return true;
    } );

    // TODO: test phi3s, maybe in three topologies for ease: 1x1x2, 1x2x1, 2x1x1
}
