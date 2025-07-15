#include <catch2/catch_test_macros.hpp>
#include <CustomEigen.hpp>
#include <CommonUtils.hpp>
#include <iostream>

TEST_CASE( "Write and read Eigen matrix to/from file" )
{
    Eigen::MatrixXd original( 5, 7 );
    original.setRandom();
    const std::string filename = "test_matrix.txt";
    writeToFile( original, filename );
    Eigen::MatrixXd read = readFromFile( filename );
    REQUIRE( read.rows() == original.rows() );
    REQUIRE( read.cols() == original.cols() );
    std::cout << "Original matrix:\n" << original << "\n";
    std::cout << "Read matrix:\n" << read << "\n";
    REQUIRE( util::equals( original, read, 1e-10 ) );
}