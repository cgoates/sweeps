#include <ParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>

namespace param
{
    BaryCoordIsZeroVec parentDomainBoundary( const ParametricAtlas& atlas, const topology::Cell& cell )
    {
        if( cell.dim() == 0 ) return atlas.parentPoint( cell ).mBaryCoordIsZero;
        else
        {
            BaryCoordIsZeroVec out( numTotalCoordinates( atlas.parentDomain( topology::Cell( cell.dart(), atlas.cmap().dim() ) ) ), true );
            iterateAdjacentCellsOfRestrictedCell( atlas.cmap(), cell, atlas.cmap().dim(), 0, [&]( const topology::Vertex& v ) {
                const ParentPoint ppt = atlas.parentPoint( v );
                out = join( out, ppt.mBaryCoordIsZero );
                return true;
            } );
            return out;
        }
    }
}