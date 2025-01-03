#pragma once
#include <CombinatorialMap.hpp>

namespace reparam
{
    int lowerLinkEulerCharacteristic( const topology::CombinatorialMap& map,
                                      const topology::Vertex& v,
                                      const std::function<double( const topology::Vertex& )>& f );

    bool isCriticalPoint( const topology::CombinatorialMap& map,
                          const topology::Vertex& v,
                          const std::function<double( const topology::Vertex& )>& f );
} // namespace reparam