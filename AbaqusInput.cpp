#include <AbaqusInput.hpp>
#include <fstream>
#include <map>
#include <Logging.hpp>
#include <algorithm>
#include <sstream>

#define LOG_INP_LOADING 0

namespace io
{

    SweepInput loadINPFile( const std::string& filename,
                            const std::string& zero_bc_label,
                            const std::string& one_bc_label )
    {
        std::ifstream read( filename );

        enum class ReadState
        {
            Initial,
            Vertices,
            ZeroBCs,
            OneBCs,
            Tets
        };

        std::vector<Eigen::Vector3d> points;
        std::map<size_t, size_t> points_reindexing;
        std::vector<Simplex> tets;
        std::set<VertexId> zero_bcs;
        std::set<VertexId> one_bcs;

        ReadState state = ReadState::Initial;
        for( std::string line; std::getline( read, line ); )
        {
            if( line.rfind( "*", 0 ) == 0 )
            {
                if( line.rfind( "*NODE", 0 ) == 0 )
                {
                    state = ReadState::Vertices;
                }
                else if( line.rfind( "*ELEMENT", 0 ) == 0 )
                {
                    if( line.rfind( " type=C3D4," ) != line.npos )
                        state = ReadState::Tets;
                    else if( line.rfind( " type=CPS3," ) != line.npos )
                    {
                        if( line.rfind( "ELSET=" + zero_bc_label ) == line.size() - zero_bc_label.size() - 6 )
                            state = ReadState::ZeroBCs;
                        else if( line.rfind( "ELSET=" + one_bc_label ) == line.size() - one_bc_label.size() - 6 )
                            state = ReadState::OneBCs;
                        else
                            state = ReadState::Initial;
                    }
                    else
                        state = ReadState::Initial;
                }
                else if( line.rfind( "*ELSET", 0 ) == 0 )
                {
                    state = ReadState::Initial;
                }
                // Otherwise it's just a comment
            }
            else if( state == ReadState::Vertices )
            {
                std::replace_if(
                    std::begin( line ), std::end( line ), []( char x ) { return x == ','; }, ' ' );

                size_t vert_id;
                double x, y, z;
                std::stringstream ss( line );
                ss >> vert_id >> x >> y >> z;

                const Eigen::Vector3d point( x, y, z );
                points_reindexing.emplace( vert_id, points.size() );
                points.push_back( point );
                LOG( LOG_INP_LOADING ) << point.transpose() << std::endl;
            }
            else if( state == ReadState::Tets )
            {
                std::replace_if(
                    std::begin( line ), std::end( line ), []( char x ) { return x == ','; }, ' ' );

                size_t tet_id, v0, v1, v2, v3;
                std::stringstream ss( line );
                ss >> tet_id >> v0 >> v1 >> v2 >> v3;

                tets.emplace_back( points_reindexing.at( v0 ),
                                   points_reindexing.at( v1 ),
                                   points_reindexing.at( v2 ),
                                   points_reindexing.at( v3 ) );
                LOG( LOG_INP_LOADING ) << tets.back() << std::endl;
            }
            else if( state != ReadState::Initial )
            {
                std::replace_if(
                    std::begin( line ), std::end( line ), []( char x ) { return x == ','; }, ' ' );
                size_t tri_id, v0, v1, v2;
                std::stringstream ss( line );
                ss >> tri_id >> v0 >> v1 >> v2;
                auto& bc_set = state == ReadState::ZeroBCs ? zero_bcs : one_bcs;
                bc_set.emplace( points_reindexing.at( v0 ) );
                bc_set.emplace( points_reindexing.at( v1 ) );
                bc_set.emplace( points_reindexing.at( v2 ) );
            }
        }

        return SweepInput::fromSets( { tets, points }, zero_bcs, one_bcs );
    }

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

    std::vector<BarycentricPoint> loadBaryCoords( const std::string& bary_filename_in, const std::set<size_t>& edge_ii )
    {
        std::vector<BarycentricPoint> output;

        const auto process_vertex_line = [&output]( const std::string& line ) {
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
        size_t n_edges = 0;
        for( std::string line; std::getline( read, line ); )
        {
            if( std::find( edge_ii.begin(), edge_ii.end(), n_edges ) != edge_ii.end() and line.rfind( "End Location", 0 ) == 0 )
            {
                std::getline( read, end_location_line );
            }
            if( std::find( edge_ii.begin(), edge_ii.end(), n_edges ) != edge_ii.end() and line.rfind( "v ", 0 ) == 0 )
            {
                process_vertex_line( line );
            }
            else if( line.rfind( "#contouredge", 0 ) == 0 )
            {
                if( end_location_line.size() > 0 )
                {
                    process_vertex_line( end_location_line );
                    end_location_line = "";
                    // FIXME: If we load all lines, the last end location will be dropped...
                }
                n_edges++;
                // if( n_edges > max_edges ) break;
            }
        }

        return output;
    }
} // namespace io