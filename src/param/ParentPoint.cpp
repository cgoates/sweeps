#include <ParentDomain.hpp>
#include <ParentPoint.hpp>
#include <CommonUtils.hpp>

namespace param
{
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
        return expandedCoordinates( pt.mDomain, pt.mPoint );
    }

    ParentPoint average( const ParentPoint& pt1, const ParentPoint& pt2 )
    {
        if( pt1.mDomain != pt2.mDomain ) throw std::runtime_error( "Cannot average two parent points from different domains" );
        return ParentPoint( pt1.mDomain, 0.5 * ( pt1.mPoint + pt2.mPoint ), join( pt1.mBaryCoordIsZero, pt2.mBaryCoordIsZero ) );
    }

    ParentPoint tensorProduct( const ParentPoint& pt1, const ParentPoint& pt2 )
    {
        BaryCoordIsZeroVec zero_vec = pt1.mBaryCoordIsZero;
        for( const bool b : pt2.mBaryCoordIsZero ) zero_vec.push_back( b );
        return ParentPoint(
            tensorProduct( pt1.mDomain, pt2.mDomain ),
            ( Vector3dMax( pt1.mPoint.size() + pt2.mPoint.size() ) << pt1.mPoint, pt2.mPoint ).finished(),
            zero_vec );
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