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

    void KnotVector::iterateUniqueKnots( const std::function<void( const double, const size_t )>& iter ) const
    {
        for( const auto& pr : mKnots ) iter( pr.first, pr.second );
    }

    void KnotVector::duplicate( const size_t knot_ii )
    {
        size_t knot_accumulator = 0;
        for( auto& pr : mKnots )
        {
            knot_accumulator += pr.second;
            if( knot_accumulator > knot_ii )
            {
                pr.second++;
                return;
            }
        }
        throw std::runtime_error( "knot_ii provided is greater than the size of the knot vector" );
    }

    Eigen::VectorXd parametricLengths( const KnotVector& kv )
    {
        std::vector<double> unique_knots = kv.uniqueKnots();
        if( unique_knots.size() < 2 ) return Eigen::VectorXd();
        return Eigen::Map<Eigen::VectorXd>( std::next( unique_knots.data() ), unique_knots.size() - 1 ) -
               Eigen::Map<Eigen::VectorXd>( unique_knots.data(), unique_knots.size() - 1 );
    }

    using SparseMatrixXd = Eigen::SparseMatrix<double>;
    SparseMatrixXd globalExtractionOp( const KnotVector& kv, const size_t degree )
    {
        const size_t num_funcs = kv.size() - degree - 1;
        SparseMatrixXd C( num_funcs, num_funcs );
        C.setIdentity();
        KnotVector knots = kv;

        size_t r_ii = 0;
        const auto add_knot = [&]( const double knot_to_add, const size_t interval_ii ) {
            SparseMatrixXd C_ii( num_funcs + r_ii, num_funcs + r_ii + 1 );
            C_ii.reserve( ( Eigen::VectorXi( num_funcs + r_ii + 1 ) << 1, Eigen::VectorXi( num_funcs + r_ii - 1, 2 ), 1 ).finished() );
            for( size_t i = 0; i < num_funcs + r_ii + 1; i++ )
            {
                const double alpha = [&]() {
                    if( i <= interval_ii - degree )
                        return 1.0;
                    else if( i < interval_ii )// Would need to be <= if we weren't adding an existing knot
                        return ( knot_to_add - knots.knot( i ) ) / ( knots.knot( i + degree ) - knots.knot( i ) );
                    else
                        return 0.0;
                }();

                if( i > interval_ii - degree )
                    C_ii.insert( i - 1, i ) = 1 - alpha;
                if( i < num_funcs + r_ii and i < interval_ii )
                    C_ii.insert( i, i ) = alpha;
            }
            C = ( C * C_ii ).eval();
            knots.duplicate( interval_ii );
            r_ii++;
        };

        size_t knot_accumulator = 0;
        kv.iterateUniqueKnots( [&]( const double knot_val, size_t mult ) {
            knot_accumulator += mult;
            while( mult < degree )
            {
                add_knot( knot_val, knot_accumulator - 1 );
                mult++;
                knot_accumulator++;
            }
        } );

        C.makeCompressed();
        return C;
    }
}