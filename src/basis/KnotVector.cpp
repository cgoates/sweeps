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

    KnotVector::KnotVector( const std::vector<std::pair<double, size_t>>& knot_multiplicities )
        : mKnots( knot_multiplicities )
    {}

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
        const std::vector<double> unique_knots = kv.uniqueKnots();
        if( unique_knots.size() < 2 ) return Eigen::VectorXd();
        return Eigen::Map<const Eigen::VectorXd>( std::next( unique_knots.data() ), unique_knots.size() - 1 ) -
               Eigen::Map<const Eigen::VectorXd>( unique_knots.data(), unique_knots.size() - 1 );
    }

    size_t numElements( const KnotVector& kv )
    {
        return kv.uniqueKnotMultiplicities().size() - 1;
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
            C_ii.reserve( ( Eigen::VectorXi( num_funcs + r_ii + 1 ) << 1, Eigen::VectorXi::Constant( num_funcs + r_ii - 1, 2 ), 1 ).finished() );
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

    KnotVector reducedOrder( const KnotVector& kv )
    {
        std::vector<std::pair<double, size_t>> knots = kv.uniqueKnotMultiplicities();

        if( knots.front().second < 2 or knots.back().second < 2 )
            throw std::runtime_error( "Reducing order of piecewise constant, floating, or periodic splines is unsupported" );

        knots.front().second--;
        knots.back().second--;

        return KnotVector( knots );
    }

    Eigen::VectorXd grevillePoints( const KnotVector& kv, const size_t degree )
    {
        const auto average_p_points = [&]( const size_t i ) {
            double ave = 0;
            for( size_t j = 0; j < degree; j++ ) ave += kv.knot( i + j ) / degree;
            return ave;
        };

        Eigen::VectorXd out( kv.size() - degree - 1 );
        for( size_t i = 1, e = kv.size() - degree; i < e; i++ )
        {
            out( i - 1 ) = average_p_points( i );
        }
        return out;
    }
}