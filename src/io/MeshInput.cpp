#include <MeshInput.hpp>
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
        if( read.fail() )
        {
            throw std::runtime_error( "Could not read file" );
        }

        enum class ReadState
        {
            Initial,
            Vertices,
            ZeroBCs,
            OneBCs,
            ZeroBCsNSet,
            OneBCsNSet,
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
                    if( line.rfind( " type=C3D4," ) != line.npos or
                        line.rfind( " TYPE=C3D4," ) != line.npos )
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
                else if( line.rfind( "*NSET", 0 ) == 0 )
                {
                    if( line.rfind( "NSET=" + zero_bc_label ) != line.npos )
                    {
                        state = ReadState::ZeroBCsNSet;
                    }
                    else if( line.rfind( "NSET=" + one_bc_label ) != line.npos )
                    {
                        state = ReadState::OneBCsNSet;
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
            else if( state == ReadState::ZeroBCs or state == ReadState::OneBCs )
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
            else if( state == ReadState::ZeroBCsNSet or state == ReadState::OneBCsNSet )
            {
                std::replace_if(
                    std::begin( line ), std::end( line ), []( char x ) { return x == ','; }, ' ' );
                size_t v;
                std::stringstream ss( line );
                auto& bc_set = state == ReadState::ZeroBCsNSet ? zero_bcs : one_bcs;
                while( ss >> v )
                {
                    bc_set.emplace( points_reindexing.at( v ) );
                }
            }
        }

        return SweepInput::fromSets( { tets, points }, zero_bcs, one_bcs );
    }

    std::pair<std::vector<SmallVector<VertexId, 4>>, std::vector<Eigen::Vector3d>>
        loadOBJFile( const std::string& filename )
    {
        std::ifstream read( filename );
        if( read.fail() )
        {
            throw std::runtime_error( "Could not read file" );
        }

        std::vector<Eigen::Vector3d> points;
        std::vector<SmallVector<VertexId, 4>> faces;

        for( std::string line; std::getline( read, line ); )
        {
            // Skip empty lines and comments
            if( line.empty() || line[0] == '#' )
            {
                continue;
            }

            std::stringstream ss( line );
            std::string prefix;
            ss >> prefix;

            if( prefix == "v" )
            {
                // Vertex line: v x y z
                double x, y, z;
                ss >> x >> y >> z;

                points.emplace_back( x, y, z );
            }
            else if( prefix == "f" )
            {
                // Face line: f v1 v2 v3 [v4] (possibly with texture/normal indices)
                // Parse vertex indices, handling potential texture/normal coordinates (v/vt/vn format)
                auto parseVertexIndex = []( const std::string& vertex_str ) -> size_t {
                    size_t slash_pos = vertex_str.find( '/' );
                    std::string index_str = vertex_str.substr( 0, slash_pos );
                    return std::stoul( index_str ) - 1; // OBJ uses 1-based indexing
                };

                SmallVector<VertexId, 4> face;
                std::string vertex_str;

                // Read all vertices for this face
                while( ss >> vertex_str )
                {
                    if( face.size() == 4 )
                    {
                        throw std::runtime_error( "Face has more than 4 vertices, which is not supported" );
                    }

                    size_t vertex_idx = parseVertexIndex( vertex_str );

                    if( vertex_idx >= points.size() )
                    {
                        throw std::runtime_error( "Invalid vertex index in face definition" );
                    }

                    face.push_back( vertex_idx );
                }

                if( face.size() < 3 )
                {
                    throw std::runtime_error( "Face must have at least 3 vertices" );
                }

                faces.push_back( face );
            }
        }

        return std::make_pair( faces, points );
    }

} // namespace io