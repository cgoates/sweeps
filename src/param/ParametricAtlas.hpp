#pragma once
#include <CombinatorialMap.hpp>
#include <ParentDomain.hpp>

namespace param
{
    class ParametricAtlas
    {
        virtual const topology::CombinatorialMap& cmap() const = 0;
        virtual const ParentDomain& parentDomain( const topology::Cell& c ) const = 0;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const = 0;
    };
}