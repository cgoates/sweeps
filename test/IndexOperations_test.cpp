#include <catch2/catch_test_macros.hpp>
#include <IndexOperations.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

using namespace util;
TEST_CASE( "iterateTensorProduct with permutations" )
{
    const auto test_iterate_tp = []( const IndexVec& lengths, const IndexVec& order, const std::vector<IndexVec>& expected ) {
        std::vector<IndexVec> actual;
        actual.reserve( expected.size() );

        iterateTensorProduct( lengths, order, [&]( const auto& iv ) {
            actual.push_back( iv );
        } );

        CHECK( expected == actual );
    };

    test_iterate_tp( { 4, 3 },
                     { 0, 1 },
                     { { 0, 0 },
                       { 1, 0 },
                       { 2, 0 },
                       { 3, 0 },
                       { 0, 1 },
                       { 1, 1 },
                       { 2, 1 },
                       { 3, 1 },
                       { 0, 2 },
                       { 1, 2 },
                       { 2, 2 },
                       { 3, 2 } } );

    test_iterate_tp( { 4, 3 },
                     { 1, 0 },
                     { { 0, 0 },
                       { 0, 1 },
                       { 0, 2 },
                       { 1, 0 },
                       { 1, 1 },
                       { 1, 2 },
                       { 2, 0 },
                       { 2, 1 },
                       { 2, 2 },
                       { 3, 0 },
                       { 3, 1 },
                       { 3, 2 } } );

    test_iterate_tp( { 4, 2, 3 }, { 0, 1, 2 }, { { 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 }, { 3, 0, 0 }, { 0, 1, 0 },
                                                 { 1, 1, 0 }, { 2, 1, 0 }, { 3, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 },
                                                 { 2, 0, 1 }, { 3, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 }, { 2, 1, 1 },
                                                 { 3, 1, 1 }, { 0, 0, 2 }, { 1, 0, 2 }, { 2, 0, 2 }, { 3, 0, 2 },
                                                 { 0, 1, 2 }, { 1, 1, 2 }, { 2, 1, 2 }, { 3, 1, 2 } } );

    test_iterate_tp( { 4, 2, 3 }, { 0, 2, 1 }, { { 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 }, { 3, 0, 0 }, { 0, 0, 1 },
                                                 { 1, 0, 1 }, { 2, 0, 1 }, { 3, 0, 1 }, { 0, 0, 2 }, { 1, 0, 2 },
                                                 { 2, 0, 2 }, { 3, 0, 2 }, { 0, 1, 0 }, { 1, 1, 0 }, { 2, 1, 0 },
                                                 { 3, 1, 0 }, { 0, 1, 1 }, { 1, 1, 1 }, { 2, 1, 1 }, { 3, 1, 1 },
                                                 { 0, 1, 2 }, { 1, 1, 2 }, { 2, 1, 2 }, { 3, 1, 2 } } );

    test_iterate_tp( { 4, 2, 3 }, { 2, 0, 1 }, { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 2 }, { 1, 0, 0 }, { 1, 0, 1 },
                                                 { 1, 0, 2 }, { 2, 0, 0 }, { 2, 0, 1 }, { 2, 0, 2 }, { 3, 0, 0 },
                                                 { 3, 0, 1 }, { 3, 0, 2 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, 1, 2 },
                                                 { 1, 1, 0 }, { 1, 1, 1 }, { 1, 1, 2 }, { 2, 1, 0 }, { 2, 1, 1 },
                                                 { 2, 1, 2 }, { 3, 1, 0 }, { 3, 1, 1 }, { 3, 1, 2 } } );
}

TEST_CASE( "iterateTensorProduct with permutations and flips" )
{
    const auto test_iterate_tp = []( const IndexVec& lengths,
                                     const IndexVec& order,
                                     const SmallVector<std::variant<bool, size_t>, 3>& dir,
                                     const std::vector<IndexVec>& expected ) {
        std::vector<IndexVec> actual;
        actual.reserve( expected.size() );

        iterateTensorProduct( lengths, order, dir, [&]( const auto& iv ) {
            actual.push_back( iv );
        } );

        CHECK( expected == actual );
    };

    test_iterate_tp( { 4, 2, 3 },
                     { 1, 2, 0 },
                     { true, size_t( 6 ), false },
                     { { 0, 6, 2 },
                       { 0, 6, 1 },
                       { 0, 6, 0 },
                       { 1, 6, 2 },
                       { 1, 6, 1 },
                       { 1, 6, 0 },
                       { 2, 6, 2 },
                       { 2, 6, 1 },
                       { 2, 6, 0 },
                       { 3, 6, 2 },
                       { 3, 6, 1 },
                       { 3, 6, 0 } } );

    test_iterate_tp( { 3, 3, 3 },
                     { 1, 0, 2 },
                     { false, false, size_t( 2 ) },
                     { { 2, 2, 2 },
                       { 2, 1, 2 },
                       { 2, 0, 2 },
                       { 1, 2, 2 },
                       { 1, 1, 2 },
                       { 1, 0, 2 },
                       { 0, 2, 2 },
                       { 0, 1, 2 },
                       { 0, 0, 2 } } );
}

TEST_CASE( "iterateTensorProductSynchronized" )
{
    const std::vector<std::pair<IndexVec, IndexVec>> expected{
        { { 0, 0, 0 }, { 0, 3, 1 } },
        { { 1, 0, 0 }, { 0, 2, 1 } },
        { { 2, 0, 0 }, { 0, 1, 1 } },
        { { 3, 0, 0 }, { 0, 0, 1 } },
        { { 0, 0, 1 }, { 1, 3, 1 } },
        { { 1, 0, 1 }, { 1, 2, 1 } },
        { { 2, 0, 1 }, { 1, 1, 1 } },
        { { 3, 0, 1 }, { 1, 0, 1 } },
        { { 0, 0, 2 }, { 2, 3, 1 } },
        { { 1, 0, 2 }, { 2, 2, 1 } },
        { { 2, 0, 2 }, { 2, 1, 1 } },
        { { 3, 0, 2 }, { 2, 0, 1 } },
    };
    std::vector<std::pair<IndexVec, IndexVec>> actual;
    actual.reserve( expected.size() );

    util::iterateTensorProductSynchronized(
        { 4, 2, 3 },
        { 3, 4, 2 },
        { 0, 1, 2 },
        { 1, 2, 0 },
        { true, size_t( 0 ), true },
        { true, false, size_t{ 1 } },
        [&]( const util::IndexVec& a, const util::IndexVec& b ) { actual.push_back( { a, b } ); } );

    CHECK( expected == actual );
}