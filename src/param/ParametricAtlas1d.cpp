#include <ParametricAtlas1d.hpp>
#include <CombinatorialMapMethods.hpp>

namespace param
{
    ParametricAtlas1d::ParametricAtlas1d( const std::shared_ptr<const topology::CombinatorialMap1d>& cmap,
                                          const Eigen::VectorXd& lengths )
        : mMap( cmap ), mLengths( lengths )
    {}

    ParametricAtlas1d::ParametricAtlas1d( const std::shared_ptr<const topology::CombinatorialMap1d>& cmap )
        : ParametricAtlas1d( cmap, Eigen::VectorXd::Ones( topology::cellCount( *cmap, cmap->dim() ) ) )
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

    std::vector<std::pair<topology::Cell, param::ParentPoint>> parentPointsOfParamPoints(
        const std::vector<double>& values, const param::ParametricAtlas1d& pa, const double param_tol )
    {
        std::vector<std::pair<topology::Cell, param::ParentPoint>> out;
        out.reserve( values.size() );

        size_t value_idx = 0;
        double current_position = values.front();
        const double factor = ( values.back() - values.front() ) / pa.totalLength();

        const param::ParentDomain pd = param::cubeDomain( 1 );

        iterateCellsWhile( pa.cmap(), 1, [&]( const topology::Edge& c ) {
            const double interval_length = pa.parametricLengths( c )( 0 ) * factor;
            const double next_position = current_position + interval_length;

            // Process all values that fall within the current interval
            for( ; value_idx < values.size() and values.at( value_idx ) <= next_position; value_idx++ )
            {
                const auto [relative_pos, zerovec] = [&]() -> std::pair<double, param::BaryCoordIsZeroVec> {
                    if( values.at( value_idx ) <= current_position + param_tol )
                    {
                        return { 0.0, { false, true } };
                    }
                    if( values.at( value_idx ) >= next_position - param_tol )
                    {
                        return { 1.0, { true, false } };
                    }
                    return { ( values.at( value_idx ) - current_position ) / interval_length, { false, false } };
                }();

                out.emplace_back( c, param::ParentPoint( pd, Vector1d( relative_pos ), zerovec ) );
            }

            current_position = next_position;
            return true;
        } );

        for( ; value_idx < values.size(); value_idx++ )
        {
            if( values.at( value_idx ) <= current_position + param_tol )
                out.emplace_back( topology::Edge( pa.cmap().maxDartId() ),
                                  param::ParentPoint( pd, Vector1d( 1.0 ), { true, false } ) );
            else
                throw std::runtime_error( "Level set values outside of parametric domain" );
        }

        return out;
    }
}