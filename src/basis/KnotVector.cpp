#include <KnotVector.hpp>
#include <numeric>
#include <CommonUtils.hpp>

namespace basis
{
    KnotVector::KnotVector( const std::vector<double>& knots, const double parametric_tol )
    {
        if( knots.size() == 0 ) throw std::runtime_error( "Cannot create an empty knot vector" );
        double current_knot = knots.front();
        mKnots.push_back( { current_knot, 1 } );
        for( size_t i = 1; i < knots.size(); i++ )
        {
            if( util::equals( current_knot, knots.at( i ), parametric_tol ) )
            {
                mKnots.back().second++;
            }
            else
            {
                current_knot = knots.at( i );
                mKnots.push_back( { current_knot, 1 } );
            }
        }
    }

    size_t KnotVector::size() const
    {
        return std::transform_reduce(
            mKnots.begin(), mKnots.end(), 0, std::plus<>(), []( const auto& pr ) { return pr.second; } );
    }

    double KnotVector::knot( const size_t knot_ii ) const
    {
        size_t knot_accumulator = 0;
        for( const auto& pr : mKnots )
        {
            knot_accumulator += pr.second;
            if( knot_accumulator > knot_ii ) return pr.first;
        }
        throw std::runtime_error( "knot_ii provided is greater than the size of the knot vector" );
    }

    std::vector<double> KnotVector::uniqueKnots() const
    {
        std::vector<double> out;
        out.reserve( mKnots.size() );
        std::transform( mKnots.begin(), mKnots.end(), std::back_inserter( out ), []( const auto& pr ){ return pr.first; } );
        return out;
    }
}