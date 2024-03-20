#include <QuadLayoutInput.hpp>
#include <fstream>
#include <map>
#include <Logging.hpp>
#include <algorithm>
#include <sstream>

namespace io
{
    std::map<size_t, size_t> vertexCorrespondences( const SimplicialComplex& mesh3d, const std::string& filename )
    {
        std::ifstream read( filename );

        std::map<size_t, size_t> points_reindexing;

        size_t i = 0;
        for( std::string line; std::getline( read, line ); )
        {
            if( line.rfind( "v ", 0 ) == 0 )
            {
                Eigen::Vector3d v_pos;
                std::string dummy_string;
                std::stringstream ss( line );
                ss >> dummy_string >> v_pos(0) >> v_pos(1) >> v_pos(2);

                const auto it = std::find_if( mesh3d.points.begin(), mesh3d.points.end(), [&]( const Eigen::Vector3d& p ) {
                    return equals( p, v_pos, 1e-4 );
                } );
                if( it == mesh3d.points.end() )
                    throw( "Point doesn't exist!!" );

                points_reindexing.emplace( i++, it - mesh3d.points.begin() );
            }
        }

        return points_reindexing;
    }

    void rewriteBaryCoordFile( const SimplicialComplex& mesh3d,
                               const std::string& obj_filename,
                               const std::string& bary_filename_in,
                               const std::string& bary_filename_out )
    {
        const std::map<size_t, size_t> c = vertexCorrespondences( mesh3d, obj_filename );

        std::ifstream read( bary_filename_in );

        std::ofstream out_file;
        out_file.open( bary_filename_out );

        for( std::string line; std::getline( read, line ); )
        {
            if( line.rfind( "v ", 0 ) == 0 )
            {
                /*
                Split the string by spaces and check each resulting piece to see if it's a v.
                If it is then the next one needs to be changed.
                 */

                bool first = true;
                std::istringstream iss( line );
                std::string s;
                while ( getline( iss, s, ' ' ) ) {
                    if( not first ) out_file << " ";
                    out_file << s;
                    first = false;
                    if( s == "v" )
                    {
                        getline( iss, s, ' ' );
                        out_file << " " << c.at( stoi( s ) );
                    }
                }
                out_file << std::endl;
            }
            else
            {
                out_file << line << std::endl;
            }
        }
    }

    std::vector<std::vector<BarycentricPoint>> loadBaryCoords( const std::string& bary_filename_in )
    {
        std::vector<std::vector<BarycentricPoint>> output;
        output.push_back({});

        const auto process_vertex_line = []( const std::string& line, std::vector<BarycentricPoint>& output ) {
            /*
            Split the string by spaces and check each resulting piece to see if it's a v.
            If it is then the next one needs to be changed.
            */

            std::istringstream iss( line );
            std::string s;
            size_t v1, v2, v3;
            double b;
            iss >> s;
            iss >> v1;
            if( (iss >> s).eof() )
            {
                output.push_back( { Simplex( v1 ), Eigen::Matrix<double, 1, 1>( 1.0 ) } );
            }
            else
            {
                iss >> v2;
                iss >> s;
                if( s == "b" )
                {
                    iss >> b;
                    if( v2 < v1 )
                        output.push_back( {Simplex(v2, v1), Eigen::Vector2d( b, 1.0 - b ) } );
                    else
                        output.push_back( {Simplex(v1, v2), Eigen::Vector2d( 1.0 - b, b ) } );
                }
                else
                {
                    double b2, b3;
                    iss >> v3 >> s >> b >> s >> b2 >> s >> b3;

                    if( v1 > v2 )
                    {
                        std::swap( v1, v2 );
                        std::swap( b, b2 );
                    }
                    if( v1 > v3 )
                    {
                        std::swap( v1, v3 );
                        std::swap( b, b3 );
                    }
                    if( v2 > v3 )
                    {
                        std::swap( v2, v3 );
                        std::swap( b2, b3 );
                    }
                    output.push_back( { Simplex( v1, v2, v3 ), Eigen::Vector3d( b, b2, b3 ) } );
                }
            }
        };

        std::ifstream read( bary_filename_in );

        std::string end_location_line = "";
        for( std::string line; std::getline( read, line ); )
        {
            if( line.rfind( "End Location", 0 ) == 0 )
            {
                std::getline( read, end_location_line );
            }
            if( line.rfind( "v ", 0 ) == 0 )
            {
                process_vertex_line( line, output.back() );
            }
            else if( line.rfind( "#contouredge", 0 ) == 0 )
            {
                if( end_location_line.size() > 0 )
                {
                    process_vertex_line( end_location_line, output.back() );
                    end_location_line = "";
                    output.back().erase( std::unique( output.back().begin(), output.back().end() ), output.back().end() );
                    output.push_back({});
                }
            }
        }

        if( end_location_line.size() > 0 )
        {
            process_vertex_line( end_location_line, output.back() );
            output.back().erase( std::unique( output.back().begin(), output.back().end() ), output.back().end() );
        }

        return output;
    }
} // namespace io