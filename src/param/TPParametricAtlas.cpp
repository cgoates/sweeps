#include <TPParametricAtlas.hpp>
#include <utility>
#include <CombinatorialMapMethods.hpp>

namespace param
{
    TPParametricAtlas::TPParametricAtlas( const topology::TPCombinatorialMap& cmap,
                                          const ParametricAtlas& source_atlas,
                                          const ParametricAtlas1d& line_atlas )
        : mMap( cmap ), mSourceParam( source_atlas ), mLineParam( line_atlas )
    {
        if( &cmap.sourceCMap() != &source_atlas.cmap() or &cmap.lineCMap() != &line_atlas.cmap() )
            throw std::runtime_error( "Incompatible cmap and atlases provided to TPParametricAtlas" );
    }

    const ParentDomain TPParametricAtlas::parentDomain( const topology::Cell& c ) const
    {
        const auto [source_dart, line_dart, tp_pos] = mMap.unflatten( c.dart() );
        const ParentDomain source_domain = mSourceParam.parentDomain( topology::Cell( source_dart, mSourceParam.cmap().dim() ) );
        const ParentDomain line_domain = mLineParam.parentDomain( topology::Cell( line_dart, mLineParam.cmap().dim() ) );
        return tensorProduct( source_domain, line_domain );
    }

    ParentPoint TPParametricAtlas::parentPoint( const topology::Vertex& v ) const
    {
        // We do some finagling to get the appropriate points from the underlying atlases.

        const auto [source_dart, line_dart, tp_pos] = mMap.unflatten( v.dart() );

        const ParentPoint line_ppt = [&]() {
            if( std::to_underlying( tp_pos ) >= mMap.dim() )
                return ParentPoint( simplexDomain( 1 ), Eigen::Matrix<double, 1, 1>( 1.0 ), { true, false } );
            else
                return ParentPoint( simplexDomain( 1 ), Eigen::Matrix<double, 1, 1>( 0.0 ), { false, true } );
        }();

        const ParentPoint source_ppt = [&]() {
            if( mMap.dim() == 3 )
            {
                switch( tp_pos )
                {
                    case topology::TPCombinatorialMap::TPDartPos::DartPos0:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos2:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos3:
                        return mSourceParam.parentPoint( source_dart );
                    default:
                        return mSourceParam.parentPoint( phi( mSourceParam.cmap(), 1, source_dart ).value() );
                }
            }
            else
            {
                switch( tp_pos )
                {
                    case topology::TPCombinatorialMap::TPDartPos::DartPos0:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos3:
                        return ParentPoint( simplexDomain( 1 ), Eigen::Matrix<double, 1, 1>( 0.0 ), { false, true } );
                    default:
                        return ParentPoint( simplexDomain( 1 ), Eigen::Matrix<double, 1, 1>( 1.0 ), { true, false } );
                }
            }
        }();

        return tensorProduct( source_ppt, line_ppt );
    }

    Vector6dMax TPParametricAtlas::parametricLengths( const topology::Cell& c ) const
    {
        const auto [source_dart, line_dart, tp_pos] = mMap.unflatten( c.dart() );
        const Vector6dMax source_lengths = mSourceParam.parametricLengths( topology::Cell( source_dart, mSourceParam.cmap().dim() ) );
        const Vector6dMax line_lengths = mLineParam.parametricLengths( topology::Cell( line_dart, mLineParam.cmap().dim() ) );
        return ( Vector6dMax( source_lengths.size() + line_lengths.size() ) << source_lengths, line_lengths ).finished();
    }
}