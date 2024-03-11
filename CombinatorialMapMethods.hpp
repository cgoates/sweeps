#pragma once
#include <CombinatorialMap.hpp>
#include <set>

namespace topology
{
    std::optional<Dart> phi( const CombinatorialMap& map, const int phi_op, const Dart& d );

    std::optional<Dart> phi( const CombinatorialMap& map, const std::vector<int>& phi_ops, const Dart& d );

    bool iterateDartsOfCell( const CombinatorialMap& map,
                             const Cell& c,
                             const std::function<bool( const Dart& )>& callback );

    bool iterateAdjacentCells( const CombinatorialMap& map,
                               const Cell& c,
                               const uint cell_dim,
                               const std::function<bool( const Cell& )>& callback );

    bool iterateCellsWhile( const CombinatorialMap& map,
                            const uint cell_dim,
                            const std::function<bool( const Cell& )>& callback );

    bool iterateDartsWhile( const CombinatorialMap& map,
                            const std::function<bool( const Dart& )>& callback );

    size_t cellCount( const CombinatorialMap& map, const uint cell_dim );

    bool onBoundary( const CombinatorialMap& map, const Dart& d );

    bool boundaryAdjacent( const CombinatorialMap& map, const Cell& c );

    /// Returns one dart for every boundary connected component.
    /// Topologies homogeneous to a disk or ball will have one dart
    /// in the returned set, while a 2-sphere would have zero.
    std::set<Dart> boundaryComponentDarts( const CombinatorialMap& map );
} // namespace topology