#include <TriangleParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>

namespace param
{
    const ParentDomain TriangleParametricAtlas::parentDomain( const topology::Cell& c ) const
    {
        if( c.dim() != 2 ) throw std::runtime_error( "Invalid cell dimension" );
        return mParentDomain;
    }
    
    ParentPoint TriangleParametricAtlas::parentPoint( const topology::Vertex& v ) const
    {
        // Iterate the darts of the triangle, find the lowest number.
        // That one is the origin.  Then find the number of phi1s from that to the input.
        // Should be doable with three phi(-1)s.
        size_t num_phi1s = 0;
        topology::Dart lowest_id = v.dart();
        topology::Dart curr_d = v.dart();
        for( size_t curr_num_phi1s = 1; curr_num_phi1s < 3; curr_num_phi1s++ )
        {
            curr_d = phi( *mMap, -1, curr_d ).value();
            if( curr_d < lowest_id )
            {
                lowest_id = curr_d;
                num_phi1s = curr_num_phi1s;
            }
        }
        if( phi( *mMap, -1, curr_d ).value() != v.dart() )
        {
            throw std::runtime_error( "TriangleParametricAtlas only takes triangle faces with no hanging nodes!" );
        }

        BaryCoordIsZeroVec is_zero( 3, true );
        is_zero.at( num_phi1s ) = false;

        return pointOnBoundary( mParentDomain, is_zero );
    }

    Vector6dMax TriangleParametricAtlas::parametricLengths( const topology::Cell& c ) const
    {
        if( c.dim() != 2 ) throw std::runtime_error( "Invalid cell dimension" );
        return Eigen::Vector3d::Ones();
    }

}