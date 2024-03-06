#pragma once
#include <CombinatorialMap.hpp>

namespace topology
{
    std::optional<Dart> phi( const CombinatorialMap& map, const int phi_op, const Dart& d );

    std::optional<Dart> phi( const CombinatorialMap& map, const std::vector<int>& phi_ops, const Dart& d );

    bool iterateDartsOfCell( const CombinatorialMap& map,
                             const Cell& c,
                             const std::function<bool( const Dart& )>& callback );

    size_t cellCount( const CombinatorialMap& map, const uint cell_dim );
} // namespace topology