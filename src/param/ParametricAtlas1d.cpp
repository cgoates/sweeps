#include <ParametricAtlas1d.hpp>
#include <CombinatorialMapMethods.hpp>

namespace param
{
    ParametricAtlas1d::ParametricAtlas1d( const topology::CombinatorialMap1d& cmap, const Eigen::VectorXd& lengths )
        : mMap( cmap ), mLengths( lengths )
    {}

    ParametricAtlas1d::ParametricAtlas1d( const topology::CombinatorialMap1d& cmap )
        : ParametricAtlas1d( cmap, Eigen::VectorXd::Ones( topology::cellCount( cmap, cmap.dim() ) ) )
    {}

    const ParentDomain ParametricAtlas1d::parentDomain( const topology::Cell& c ) const
    {
        if( c.dim() != 1 ) throw std::runtime_error( "Invalid cell dimension" );
        return mParentDomain;
    }

    ParentPoint ParametricAtlas1d::parentPoint( const topology::Vertex& ) const
    {
        return ParentPoint( mParentDomain, Eigen::Matrix<double, 1, 1>( 0.0 ), { false, true } );
    }

    Vector6dMax ParametricAtlas1d::parametricLengths( const topology::Cell& c ) const
    {
        if( c.dim() != 1 ) throw std::runtime_error( "invalid cell dimension" );
        return Eigen::Matrix<double, 1, 1>( mLengths( c.dart().id() ) );
    }

}