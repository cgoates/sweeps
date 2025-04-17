#include <KnotVector.hpp>
#include <numeric>
#include <CommonUtils.hpp>
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>

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

    void KnotVector::insert( const double knot_to_add, const double param_tol )
    {
        for( size_t i = 0; i < mKnots.size(); i++ )
        {
            if( util::equals( mKnots.at( i ).first, knot_to_add, param_tol ) )
            {
                mKnots.at( i ).second++;
                break;
            }
            else if( mKnots.at( i ).first > knot_to_add )
            {
                mKnots.insert( std::next( mKnots.begin(), i ), std::pair<double, size_t>{ knot_to_add, 1 } );
                break;
            }
        }
    }

    std::ostream& operator<<( std::ostream& o, const KnotVector& kv )
    {
        o << "{";
        for( size_t i = 0; i < kv.size(); i++ )
            o << kv.knot( i ) << ", ";
        o << "}";
        return o;
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
    std::pair<KnotVector, SparseMatrixXd> addKnot( const KnotVector& knots,
                                                   const size_t degree,
                                                   const double knot_to_add,
                                                   const size_t interval_ii,
                                                   const double param_tol )
    {
        const size_t num_funcs = knots.size() - degree - 1;
        SparseMatrixXd C_ii( num_funcs, num_funcs + 1 );
        C_ii.reserve( ( Eigen::VectorXi( num_funcs + 1 ) << 1, Eigen::VectorXi::Constant( num_funcs - 1, 2 ), 1 ).finished() );
        for( size_t i = 0; i < num_funcs + 1; i++ )
        {
            const double alpha = [&]() {
                if( i <= interval_ii - degree )
                    return 1.0;
                else if( i <= interval_ii )// Would need to be <= if we weren't adding an existing knot
                    return ( knot_to_add - knots.knot( i ) ) / ( knots.knot( i + degree ) - knots.knot( i ) );
                else
                    return 0.0;
            }();

            if( i > interval_ii - degree )
                C_ii.insert( i - 1, i ) = 1 - alpha;
            if( i < num_funcs and i <= interval_ii )
                C_ii.insert( i, i ) = alpha;
        }

        KnotVector new_knots = knots;
        new_knots.insert( knot_to_add, param_tol );
        return { new_knots, C_ii };
    }

    SparseMatrixXd globalExtractionOp( const KnotVector& kv, const size_t degree )
    {
        const size_t num_funcs = kv.size() - degree - 1;
        SparseMatrixXd C( num_funcs, num_funcs );
        C.setIdentity();
        KnotVector knots = kv;

        size_t knot_accumulator = 0;
        kv.iterateUniqueKnots( [&]( const double knot_val, size_t mult ) {
            knot_accumulator += mult;
            while( mult < degree )
            {
                SparseMatrixXd C_ii;
                std::tie( knots, C_ii ) = addKnot( knots, degree, knot_val, knot_accumulator - 1, 1e-12 );
                C = ( C * C_ii ).eval();
                mult++;
                knot_accumulator++;
            }
        } );

        C.makeCompressed();
        return C;
    }

    SparseMatrixXd refinementOp( const KnotVector& kv_coarse, const KnotVector& kv_fine, const size_t degree, const double param_tol )
    {
        const size_t num_funcs = kv_coarse.size() - degree - 1;
        SparseMatrixXd C( num_funcs, num_funcs );
        C.setIdentity();
        KnotVector knots = kv_coarse;

        const auto& fine_knots = kv_fine.uniqueKnotMultiplicities();
        size_t knot_accumulator = 0;
        for( size_t i = 0; i < fine_knots.size(); i++ )
        {
            const auto& coarse_knots = knots.uniqueKnotMultiplicities();
            if( i >= coarse_knots.size() )
                throw std::invalid_argument( "Fine knot vector is not nested within the coarse knot vector" );

            const double coarse_knot = coarse_knots.at( i ).first;
            if( util::equals( coarse_knot, fine_knots.at( i ).first, param_tol ) )
            {
                knot_accumulator += coarse_knots.at( i ).second;
                while( knots.uniqueKnotMultiplicities().at( i ).second < fine_knots.at( i ).second )
                {
                    SparseMatrixXd C_ii;
                    std::tie( knots, C_ii ) = addKnot( knots, degree, coarse_knot, knot_accumulator - 1, param_tol );
                    C = ( C * C_ii ).eval();
                    knot_accumulator++;
                }
            }
            else if( coarse_knot > fine_knots.at( i ).first )
            {
                do
                {
                    SparseMatrixXd C_ii;
                    std::tie( knots, C_ii ) = addKnot( knots, degree, fine_knots.at( i ).first, knot_accumulator - 1, param_tol );
                    C = ( C * C_ii ).eval();
                    knot_accumulator++;
                } while ( knots.uniqueKnotMultiplicities().at( i ).second < fine_knots.at( i ).second );
            }
            else
                throw std::invalid_argument( "Fine knot vector is not nested within the coarse knot vector" );
        }

        C.makeCompressed();
        return C;
    }

    Eigen::SparseMatrix<double> refinementOp( const SmallVector<KnotVector, 3>& kvs_coarse,
                                              const SmallVector<KnotVector, 3>& kvs_fine,
                                              const SmallVector<size_t, 3> degrees,
                                              const double param_tol )
    {
        const auto refinement_op = [&]( const size_t i ) {
            return refinementOp( kvs_coarse.at( i ), kvs_fine.at( i ), degrees.at( i ), param_tol );
        };

        Eigen::SparseMatrix<double> op = refinement_op( 0 );
        for( size_t i = 1; i < degrees.size(); i++ )
        {
            op = Eigen::kroneckerProduct( refinement_op( i ), op ).eval();
        }

        return op;
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

    KnotVector integerKnotsWithNElems( const size_t n_elems, const size_t degree )
    {
        std::vector<std::pair<double, size_t>> knots;
        knots.reserve( n_elems + 1 );
        knots.push_back( { 0.0, degree + 1 } );
        for( size_t i = 1; i < n_elems; i++ ) knots.push_back( { (double)i, 1 } );
        knots.push_back( { n_elems, degree + 1 } );

        return KnotVector( knots );
    }

    KnotVector unitIntervalKnotVectorWithNElems( const size_t n_elems, const size_t degree )
    {
        std::vector<std::pair<double, size_t>> knots;
        knots.reserve( n_elems + 1 );
        knots.push_back( { 0.0, degree + 1 } );
        for( size_t i = 1; i < n_elems; i++ ) knots.push_back( { (double)i / (double)n_elems, 1 } );
        knots.push_back( { 1.0, degree + 1 } );

        return KnotVector( knots );
    }

    KnotVector dyadicRefine( const KnotVector& kv )
    {
        return nAdicRefine( kv, 2 );
    }

    KnotVector nAdicRefine( const KnotVector& kv, const size_t n )
    {
        const auto& original_unique_knots = kv.uniqueKnotMultiplicities();
        std::vector<std::pair<double, size_t>> new_unique_knots;
        new_unique_knots.push_back( original_unique_knots.front() );
        for( size_t i = 1; i < original_unique_knots.size(); i++ )
        {
            for( size_t n_ii = 1; n_ii < n; n_ii++ )
            {
                const double x = (double)n_ii / (double)n;
                new_unique_knots.push_back( { ( 1.0 - x ) * original_unique_knots.at( i - 1 ).first + x * original_unique_knots.at( i ).first, 1 } );
            }

            new_unique_knots.push_back( original_unique_knots.at( i ) );
        }
        return KnotVector( new_unique_knots );
    }
}