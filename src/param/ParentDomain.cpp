#include <ParentDomain.hpp>
#include <CommonUtils.hpp>

namespace param
{
    size_t numTotalCoordinates( const CoordinateSystem& cs ) { return cs.dim() + 1; }
    size_t numParametricLengths( const CoordinateSystem& cs ) { return cs.dim() * ( cs.dim() + 1 ) / 2; }

    ParentDomain simplexDomain( const size_t dim )
    {
        return ParentDomain( std::initializer_list<CoordinateSystem>{ CoordinateSystem( dim ) } );
    }
    ParentDomain cubeDomain( const size_t dim )
    {
        return ParentDomain( std::vector<CoordinateSystem>( dim, CoordinateSystem( 1 ) ) );
    }

    ParentDomain tensorProduct( const ParentDomain& pd1, const ParentDomain& pd2 )
    {
        SmallVector<CoordinateSystem, 3> new_cs_vec( pd1.coordinateGroups() );
        iterateGroups( pd2, [&]( const auto&, const auto&, const CoordinateSystem& cs ){ new_cs_vec.push_back( cs ); } );
        return ParentDomain( new_cs_vec );
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

    size_t numParametricLengths( const ParentDomain& pd )
    {
        return std::accumulate( pd.coordinateGroups().begin(),
                                pd.coordinateGroups().end(),
                                0,
                                []( size_t acc, const CoordinateSystem& cs ) { return acc + numParametricLengths( cs ); } );
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

    std::ostream& operator<<( std::ostream& o, const ParentDomain& pd )
    {
        o << "ParentDomain( ";
        iterateGroups( pd, [&]( const auto&, const auto&, const CoordinateSystem& cs ) { o << cs.dim() << ", "; } );
        o << ")";
        return o;
    }
}