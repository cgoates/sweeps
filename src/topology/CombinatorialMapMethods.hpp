#pragma once
#include <CombinatorialMap.hpp>
#include <set>
#include <concepts>

namespace topology
{
    template<typename M>
    concept DartMarker = requires( M m, const Dart& d )
    {
        { m.isMarked( d ) } -> std::same_as<bool>;
        { m.mark( d ) } -> std::same_as<void>;
    };

    std::optional<Dart> phi( const CombinatorialMap& map, const int phi_op, const Dart& d );

    std::optional<Dart> phi( const CombinatorialMap& map, const std::vector<int>& phi_ops, const Dart& d );

    bool iterateDartsOfCell( const CombinatorialMap& map,
                             const Cell& c,
                             const std::function<bool( const Dart& )>& callback );

    bool iterateDartsOfCell( const CombinatorialMap& map,
                             const Cell& c,
                             DartMarker auto& m,
                             const std::function<bool( const Dart& )>& callback );

    bool iterateDartsOfRestrictedCell( const CombinatorialMap& map,
                                       const Cell& c,
                                       const int restrict_dim,
                                       const std::function<bool( const Dart& )>& callback );

    bool iterateDartsOfRestrictedCell( const CombinatorialMap& map,
                                       const Cell& c,
                                       const int restrict_dim,
                                       DartMarker auto& m,
                                       const std::function<bool( const Dart& )>& callback );

    bool iterateAdjacentCells( const CombinatorialMap& map,
                               const Cell& c,
                               const uint cell_dim,
                               const std::function<bool( const Cell& )>& callback );

    bool iterateAdjacentCellsOfRestrictedCell( const CombinatorialMap& map,
                                               const Cell& c,
                                               const int restrict_dim,
                                               const uint cell_dim,
                                               const std::function<bool( const Cell& )>& callback );

    bool iterateCellsWhile( const CombinatorialMap& map,
                            const uint cell_dim,
                            const std::function<bool( const Cell& )>& callback );

    bool iterateDartsWhile( const CombinatorialMap& map,
                            const std::function<bool( const Dart& )>& callback );

    bool iterateReachableWhile( const CombinatorialMap& map,
                                const Dart& start_dart,
                                const std::vector<int>& phi_ops,
                                const std::function<bool( const Dart& )>& callback );

    size_t cellCount( const CombinatorialMap& map, const uint cell_dim );

    bool onBoundary( const CombinatorialMap& map, const Dart& d );

    bool boundaryAdjacent( const CombinatorialMap& map, const Cell& c );
    /// If the cell is boundary adjacent, returns a dart of the cell that is on the boundary.
    /// If the cell is not boundary adjacent, returns std::nullopt.
    std::optional<Dart> maybeBoundaryDart( const CombinatorialMap& map, const Cell& c );

    bool isExtraordinary( const CombinatorialMap& map, const Vertex& v );

    /// Returns one dart for every boundary connected component.
    /// Topologies homogeneous to a disk or ball will have one dart
    /// in the returned set, while a 2-sphere would have zero.
    std::set<Dart> boundaryComponentDarts( const CombinatorialMap& map );

    /// Returns the lowest dart id of a cell.
    ///
    Dart::IndexType lowestDartId( const CombinatorialMap& map, const Cell& c );

    void flood2d( const topology::CombinatorialMap& map,
                  const topology::Face& f,
                  const std::function<bool( const topology::Face& )>& stop_condition,
                  const std::function<void( const topology::Face& )>& mark_callback,
                  const std::function<void( const topology::Face& )>& callback );

    /// Floods the entire set of darts using phi 1 and 2 operations on two combinatorial maps simultaneously.
    ///
    bool synchronizedFlood2d( const CombinatorialMap& map1,
                              const CombinatorialMap& map2,
                              const Dart& d1,
                              const Dart& d2,
                              const std::function<bool( const Dart&, const Dart& )>& callback );

    IndexingFunc indexingOrError( const CombinatorialMap& map, const uint cell_dim );

    int eulerCharacteristic( const CombinatorialMap& map );
} // namespace topology