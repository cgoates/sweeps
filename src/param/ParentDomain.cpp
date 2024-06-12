#include <ParentDomain.hpp>
#include <CommonUtils.hpp>

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

    Vector6dMax expandedCoordinates( const ParentDomain& domain, const Vector3dMax& point )
    {
        Vector6dMax out = Vector6dMax::Zero( numTotalCoordinates( domain ) );

        iterateGroups( domain, [&]( const size_t group_expanded_start, const size_t group_explicit_start, const CoordinateSystem& cs ) {
            double residue = 0;
            for( size_t group_coord = 0; group_coord < cs.dim(); group_coord++ )
            {
                residue += point( group_explicit_start + group_coord );
                out( group_expanded_start + 1 + group_coord ) = point( group_explicit_start + group_coord );
            }
            out( group_expanded_start ) = 1.0 - residue;
        } );

        return out;
    }

    ParentPoint::ParentPoint( const ParentDomain& domain, const Vector3dMax& point, const BaryCoordIsZeroVec& zero_vec ) :
        mDomain( domain ), mPoint( point ), mBaryCoordIsZero( zero_vec )
    {}
    ParentPoint compressCoordinates( const ParentDomain& domain, const Vector6dMax& coords, const double is_zero_tol ){
        BaryCoordIsZeroVec zeros( coords.size(), false );
        for( Eigen::Index i = 0; i < coords.size(); i++ )
        {
            if( util::equals( coords( i ), 0.0, is_zero_tol ) ) zeros.at( i ) = true;
        }
        Vector3dMax point = Vector3dMax::Zero( dim( domain ) );
        iterateGroups( domain, [&]( const size_t expanded_start, const size_t explicit_start, const param::CoordinateSystem& cs ) {
            for( size_t ii = 0; ii < cs.dim(); ii++ )
            {
                point( explicit_start + ii ) = coords( expanded_start + 1 + ii );
            }
        } );
        return ParentPoint( domain, point, zeros );
    }

    Vector6dMax expandedCoordinates( const ParentPoint& pt )
    {
        return expandedCoordinates( pt.mDomain, pt.mPoint );
    }

    ParentPoint average( const ParentPoint& pt1, const ParentPoint& pt2 )
    {
        if( pt1.mDomain != pt2.mDomain ) throw std::runtime_error( "Cannot average two parent points from different domains" );
        return ParentPoint( pt1.mDomain, 0.5 * ( pt1.mPoint + pt2.mPoint ), join( pt1.mBaryCoordIsZero, pt2.mBaryCoordIsZero ) );
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