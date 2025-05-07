#pragma once
#include <Cell.hpp>
#include <functional>
#include <VertexPositionsFunc.hpp>

namespace topology
{
    class CombinatorialMap;
    std::vector<Edge> shortestPath( const CombinatorialMap& map,
                                    const VertexPositionsFunc& vert_positions,
                                    const Vertex& start_vertex,
                                    const std::function<bool( const Vertex& )>& stop_condition,
                                    const bool interior_only = false );

    std::vector<Edge> shortestPath( const CombinatorialMap& map,
                                    const std::function<double( const Edge& )>& edge_lengths,
                                    const Vertex& start_vertex,
                                    const std::function<bool( const Vertex& )>& stop_condition,
                                    const bool interior_only = false );
}