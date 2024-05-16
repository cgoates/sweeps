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

    void iterateGroups( const ParentDomain& pd,
                        const std::function<void( const size_t, const size_t, const CoordinateSystem& )>& callback )
    {
        size_t group_start_expanded_coord = 0;
        size_t group_start_explicit_coord = 0;
        for( const CoordinateSystem& cs : pd.coordinateGroups() )
        {
            callback( group_start_expanded_coord, group_start_explicit_coord, cs );

            group_start_expanded_coord += numTotalCoordinates( cs );
            group_start_explicit_coord += cs.dim();
        }
    }

    ParentPoint pointOnBoundary( const ParentDomain& domain, const BaryCoordIsZeroVec& is_zero )
    {
        Vector3dMax point = Vector3dMax::Zero( dim( domain ) );
        iterateGroups( domain, [&]( const size_t group_start_expanded_coord, const size_t group_start_explicit_coord, const CoordinateSystem& cs ) {
            const size_t group_num_expanded_coords = numTotalCoordinates( cs );
            const size_t group_num_explicit_coords = cs.dim();

            const size_t num_group_nonzero_coords = [&](){
                size_t out = 0;
                for( size_t i = 0; i < group_num_expanded_coords; i++ )
                {
                    if( not is_zero.at( group_start_expanded_coord + i ) ) out++;
                }
                return out;
            }();

            const double ave_coord = 1.0 / num_group_nonzero_coords;

            for( size_t i = 0; i < group_num_explicit_coords; i++ )
            {
                if( not is_zero.at( group_start_expanded_coord + 1 + i ) )
                    point( group_start_explicit_coord + i ) = ave_coord;
            }
        } );

        return ParentPoint{ domain, point, is_zero };
    }

    Vector6dMax expandedCoordinates( const ParentPoint& pt )
    {
        Vector6dMax out = Vector6dMax::Zero( numTotalCoordinates( pt.mDomain ) );

        iterateGroups( pt.mDomain, [&]( const size_t group_expanded_start, const size_t group_explicit_start, const CoordinateSystem& cs ) {
            double residue = 0;
            for( size_t group_coord = 0; group_coord < cs.dim(); group_coord++ )
            {
                residue += pt.mPoint( group_explicit_start + group_coord );
                out( group_expanded_start + 1 + group_coord ) = pt.mPoint( group_explicit_start + group_coord );
            }
            out( group_expanded_start ) = 1.0 - residue;
        } );

        return out;
    }

    BaryCoordIsZeroVec join( const BaryCoordIsZeroVec& v1, const BaryCoordIsZeroVec& v2 )
    {
        if( v1.size() != v2.size() ) throw std::runtime_error( "Cannot join two BaryCoordIsZeroVecs of different sizes" );
        BaryCoordIsZeroVec out;
        for( size_t i = 0; i < v1.size(); i++ )
        {
            out.push_back( v1.at( i ) and v2.at( i ) );
        }
        return out;
    }
}