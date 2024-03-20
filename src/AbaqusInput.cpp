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

} // namespace io