#include <Foliation.hpp>
#include <Tracing.hpp>
#include <Cell.cpp>
#include <SimplexUtilities.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <LevelSetCMap.hpp>
#include <Logging.hpp>

namespace reparam
{
    std::vector<TraceLevelSetIntersection> levelSetIntersections( const Trace& trace,
                                                                  const topology::CombinatorialMap& traced_map,
                                                                  const std::vector<double>& level_set_values )
    {
        std::vector<TraceLevelSetIntersection> out;
        out.reserve( level_set_values.size() );

        const auto& vert_values = trace.mHarmonicFuncValues;

        size_t trace_vert_ii = 1;
        for( const double value : level_set_values )
        {
            while( vert_values.at( trace_vert_ii ) < value ) trace_vert_ii++;

            const auto [s, t] =
                inverseLinear( vert_values.at( trace_vert_ii - 1 ), vert_values.at( trace_vert_ii ), value );

            const Eigen::Vector3d pt =
                s * trace.mComplex.points.at( trace_vert_ii - 1 ) + t * trace.mComplex.points.at( trace_vert_ii );

            const topology::Cell& base_cell = trace.mBaseCells.at( trace_vert_ii - 1 );

            out.push_back( { pt, indexingOrError( traced_map, base_cell.dim() )( base_cell ) } );
        }

        return out;
    }

    std::map<topology::Vertex, double>
        thetaValues( const topology::CombinatorialMap& level_set,
                     const VertexPositionsFunc& level_set_positions,
                     const std::function<size_t( const topology::Edge& )>& underlying_face_id_of_edge,
                     const TraceLevelSetIntersection& intersection )
    {
        const topology::CombinatorialMapBoundary bdry( level_set );
        const auto bdry_positions = boundaryVertexPositions( bdry, level_set_positions );

        const topology::Dart start_d = [&bdry]() {
            topology::Dart d;
            iterateDartsWhile( bdry, [&]( const topology::Dart& a ) {
                d = a;
                return false;
            } );
            return d;
        }();

        const auto is_intersection_edge = [&]( const topology::Edge& bdry_e ) {
            const topology::Edge level_e = bdry.toUnderlyingCell( bdry_e );
            return ( underlying_face_id_of_edge( level_e ) == intersection.second );
        };

        std::map<topology::Vertex, double> out;

        topology::Dart d = start_d;
        double cumulative_length = 0.0;
        std::optional<double> intersection_point = std::nullopt;
        do
        {
            out.insert( { bdry.toUnderlyingCell( topology::Vertex( d ) ), cumulative_length } );

            if( not intersection_point and is_intersection_edge( d ) )
            {
                intersection_point.emplace( cumulative_length + ( intersection.first - bdry_positions( d ) ).norm() );
            }

            cumulative_length += edgeLength( bdry, bdry_positions, d );
            d = phi( bdry, 1, d ).value();
        } while( d != start_d );

        if( not intersection_point ) throw std::runtime_error( "Trace doesn't interact with the level set!" );

        const double factor = 2 * std::numbers::pi / cumulative_length;
        for( auto& pr : out )
        {
            pr.second = ( pr.second - intersection_point.value() ) * factor;
        }

        return out;
    }
} // namespace reparam