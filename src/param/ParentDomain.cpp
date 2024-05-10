#include <ParentDomain.hpp>

namespace param
{
    size_t numTotalCoordinates( const CoordinateSystem& cs ) { return cs.dim() + 1; }

    ParentDomain simplexDomain( const size_t dim )
    {
        return ParentDomain( std::initializer_list<CoordinateSystem>{ CoordinateSystem( dim ) } );
    }
    ParentDomain cubeDomain( const size_t dim )
    {
        return ParentDomain( std::vector<CoordinateSystem>( dim, CoordinateSystem( 1 ) ) );
    }

    size_t numGroups( const ParentDomain& pd )
    {
        return pd.coordinateGroups().size();
    }

    size_t dim( const ParentDomain& pd )
    {
        return std::accumulate( pd.coordinateGroups().begin(),
                                pd.coordinateGroups().end(),
                                0,
                                []( size_t acc, const CoordinateSystem& cs ) { return acc + cs.dim(); } );
    }

    size_t numTotalCoordinates( const ParentDomain& pd )
    {
        return std::accumulate( pd.coordinateGroups().begin(),
                                pd.coordinateGroups().end(),
                                0,
                                []( size_t acc, const CoordinateSystem& cs ) { return acc + numTotalCoordinates( cs ); } );
    }
}