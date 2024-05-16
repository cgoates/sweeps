#pragma once
#include <CombinatorialMap.hpp>
#include <ParentDomain.hpp>

namespace param
{
    class ParametricAtlas
    {
        public:
        virtual const topology::CombinatorialMap& cmap() const = 0;
        virtual const ParentDomain& parentDomain( const topology::Cell& c ) const = 0;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const = 0;
    };

    BaryCoordIsZeroVec parentDomainBoundary( const ParametricAtlas& atlas, const topology::Cell& cell );
}