#pragma once
#include <CombinatorialMap.hpp>
#include <ParentPoint.hpp>

namespace param
{
    class ParametricAtlas
    {
        public:
        virtual const topology::CombinatorialMap& cmap() const = 0;
        virtual const ParentDomain parentDomain( const topology::Cell& c ) const = 0;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const = 0;
        virtual Vector6dMax parametricLengths( const topology::Cell& c ) const = 0;
    };

    BaryCoordIsZeroVec parentDomainBoundary( const ParametricAtlas& atlas, const topology::Cell& cell );

    size_t parametricLengthIndexAlongEdge( const ParametricAtlas& atlas, const topology::Edge& e );

    topology::Vertex originVertex( const ParametricAtlas& atlas, const topology::Cell& c );

    SmallVector<std::pair<size_t, bool>, 3> coordinateTransform( const ParametricAtlas& atlas, const topology::Cell& c );

    /// A coordinate frame is a triple of parent points corresponding to the (0,0), (1,0), and (0,1) corners of a face,
    /// if the input dart is taken as aligned with the s+ direction.  If reverse_dart is true, then the dart is taken as
    /// pointing in the s- direction instead.
    SmallVector<ParentPoint, 4> getFrame( const ParametricAtlas& atlas, const topology::Cell& f, const bool reverse_dart = false );
}