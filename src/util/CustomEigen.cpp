#include <CustomEigen.hpp>
#include <fstream>

void writeToFile( const Eigen::MatrixXd& data, const std::string& filename, const int precision )
{
    std::ofstream out( filename );
    if( !out.is_open() ) return;

    out << data.rows() << " " << data.cols() << "\n";
    out << std::scientific;
    out.precision( precision );

    for( int i = 0; i < data.rows(); ++i )
    {
        for( int j = 0; j < data.cols(); ++j )
        {
            out << data( i, j );
            if( j + 1 < data.cols() ) out << " ";
        }
        out << "\n";
    }

    out.close();
}

Eigen::MatrixXd readFromFile( const std::string& filename )
{
    std::ifstream in( filename );
    if( !in.is_open() ) return Eigen::MatrixXd();

    int rows, cols;
    in >> rows >> cols;

    Eigen::MatrixXd data( rows, cols );
    for( int i = 0; i < rows; ++i )
    {
        for( int j = 0; j < cols; ++j )
        {
            in >> data( i, j );
        }
    }

    in.close();
    return data;
}