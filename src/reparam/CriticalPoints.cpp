#include <CriticalPoints.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>

using namespace topology;

static constexpr bool LOG_LOWER_LINK_EULER = false;

namespace reparam
{
    int lowerLinkEulerCharacteristic( const topology::CombinatorialMap& map,
                                      const topology::Vertex& v,
                                      const std::function<double( const topology::Vertex& )>& f )
    {
        int faces = 0;
        int edges = 0;
        int vertices = 0;

        const double val = f( v );

        LocalCellMarker v_mark( 0 );
        LocalCellMarker e_mark( 1 );
        const auto update_count = [&]( const Face& face ) {
            bool all_lower = true;
            iterateDartsOfCell( map, face, [&]( const Dart& d ) {
                if( f( d ) >= val ) all_lower = false;
                else
                {
                    if( not v_mark.isMarked( Vertex( d ) ) )
                    {
                        v_mark.mark( map, Vertex( d ) );
                        vertices++;
                    }
                    if( f( phi( map, 1, d ).value() ) < val and
                        not e_mark.isMarked( Edge( d ) ) )
                    {
                        e_mark.mark( map, Edge( d ) );
                        edges++;
                    }
                }
                return true;
            } );
            if( all_lower ) faces++;
        };

        iterateAdjacentCells( map, v, 3, [&]( const topology::Volume& v ) {
            const Face link_f( phi( map, { 1, 2 }, v.dart() ).value() );
            update_count( link_f );
            return true;
        } );

        LOG( LOG_LOWER_LINK_EULER ) << "Vertices: " << vertices << " Edges: " << edges << " Faces: " << faces << std::endl;

        return vertices - edges + faces;
    }

    bool isCriticalPoint( const topology::CombinatorialMap& map,
                          const topology::Vertex& v,
                          const std::function<double( const topology::Vertex& )>& f )
    {
        return lowerLinkEulerCharacteristic( map, v, f ) != 1;
    }
}