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

    ParentPoint pointOnBoundary( const ParentDomain& domain, const BaryCoordIsZeroVec& is_zero )
    {
        Vector3dMax point = Vector3dMax::Zero( dim( domain ) );
        size_t group_start_expanded_coord = 0;
        size_t group_start_explicit_coord = 0;
        for( const CoordinateSystem& cs : domain.coordinateGroups() )
        {
            const size_t group_num_expanded_coords = numTotalCoordinates( cs );
            const size_t group_num_explicit_coords = cs.dim();

            const size_t num_group_nonzero_coords =
                std::accumulate( std::next( is_zero.begin(), group_start_expanded_coord ),
                                 std::next( is_zero.begin(), group_start_expanded_coord + group_num_expanded_coords ),
                                 0,
                                 [&]( const bool& zero_coord, const size_t& n_nonzeros ) {
                                     return n_nonzeros + ( zero_coord ? 0 : 1 );
                                 } );

            const double ave_coord = 1.0 / num_group_nonzero_coords;

            for( size_t i = 0; i < group_num_explicit_coords; i++ )
            {
                point( group_start_explicit_coord + 1 ) = ave_coord;
            }

            group_start_expanded_coord += group_num_expanded_coords;
            group_start_explicit_coord += group_num_explicit_coords;
        }

        return ParentPoint{ domain, point, is_zero };
    }
}