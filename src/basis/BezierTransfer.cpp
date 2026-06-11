#include <BezierTransfer.hpp>
#include <CommonUtils.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace basis
{
    namespace
    {
        long double logBinomial( const size_t n, const size_t k )
        {
            if( k > n ) return -std::numeric_limits<long double>::infinity();

            return std::lgamma( static_cast<long double>( n ) + 1.0L ) -
                   std::lgamma( static_cast<long double>( k ) + 1.0L ) -
                   std::lgamma( static_cast<long double>( n - k ) + 1.0L );
        }

        double duplicateTolerance( const double param_tol )
        {
            return std::max( 100.0 * std::numeric_limits<double>::epsilon(), param_tol );
        }
    }

    std::vector<LocalExtraction> localExtractions( const KnotVector& kv, const size_t degree )
    {
        using SparseMatrixXd = Eigen::SparseMatrix<double>;

        const SparseMatrixXd C = globalExtractionOp( kv, degree );
        const std::vector<size_t> column_offsets = elementBezierColumnOffsets( kv, degree );
        std::vector<LocalExtraction> out;
        out.reserve( numElements( kv ) );

        for( size_t elem_ii = 0, n_elems = numElements( kv ); elem_ii < n_elems; elem_ii++ )
        {
            if( column_offsets.at( elem_ii ) + degree >= static_cast<size_t>( C.cols() ) )
                throw std::runtime_error( "Element Bezier column offset exceeds global extraction operator size." );

            std::set<FunctionId> unique_rows;
            for( size_t col_ii = 0; col_ii <= degree; col_ii++ )
            {
                for( SparseMatrixXd::InnerIterator it( C, col_ii + column_offsets.at( elem_ii ) ); it; ++it )
                {
                    unique_rows.insert( FunctionId( it.row() ) );
                }
            }

            LocalExtraction local;
            local.connectivity = std::vector<FunctionId>( unique_rows.begin(), unique_rows.end() );
            local.extraction = Eigen::MatrixXd::Zero( local.connectivity.size(), degree + 1 );

            for( size_t col = 0; col <= degree; col++ )
            {
                for( SparseMatrixXd::InnerIterator it( C, col + column_offsets.at( elem_ii ) ); it; ++it )
                {
                    const auto row_it =
                        std::find( local.connectivity.begin(), local.connectivity.end(), FunctionId( it.row() ) );
                    if( row_it == local.connectivity.end() )
                        throw std::runtime_error( "Local extraction connectivity is missing an extraction row." );

                    const Eigen::Index row = std::distance( local.connectivity.begin(), row_it );
                    local.extraction( row, col ) = it.value();
                }
            }

            out.push_back( local );
        }

        return out;
    }

    std::vector<size_t> elementBezierColumnOffsets( const KnotVector& kv, const size_t degree )
    {
        const std::vector<std::pair<double, size_t>> knots = kv.uniqueKnotMultiplicities();
        std::vector<size_t> offsets;
        offsets.reserve( knots.size() > 0 ? knots.size() - 1 : 0 );
        offsets.push_back( 0 );

        for( size_t i = 1; i + 1 < knots.size(); i++ )
        {
            if( knots.at( i ).second > degree + 1 )
                throw std::runtime_error( "Interior knot multiplicity exceeds degree + 1." );

            offsets.push_back( offsets.back() + std::max( knots.at( i ).second, degree ) );
        }

        return offsets;
    }

    Eigen::MatrixXd bernsteinDegreeElevationMatrix( const size_t source_degree, const size_t target_degree )
    {
        if( target_degree < source_degree ) throw std::invalid_argument( "Cannot degree elevate to a lower degree" );
        if( target_degree == source_degree ) return Eigen::MatrixXd::Identity( source_degree + 1, source_degree + 1 );

        const size_t elevation_amount = target_degree - source_degree;
        Eigen::MatrixXd out = Eigen::MatrixXd::Zero( source_degree + 1, target_degree + 1 );

        for( size_t i = 0; i <= source_degree; i++ )
        {
            for( size_t j = i; j <= i + elevation_amount; j++ )
            {
                const long double log_value =
                    logBinomial( source_degree, i ) +
                    logBinomial( elevation_amount, j - i ) -
                    logBinomial( target_degree, j );
                out( i, j ) = static_cast<double>( std::exp( log_value ) );
            }
        }

        return out;
    }

    KnotVector degreeElevatedKnotVector( const KnotVector& kv, const size_t elevation_amount )
    {
        std::vector<std::pair<double, size_t>> knots = kv.uniqueKnotMultiplicities();
        for( auto& knot_mult : knots )
        {
            knot_mult.second += elevation_amount;
        }

        return KnotVector( knots );
    }

    Eigen::SparseMatrix<double> degreeElevationOp( const KnotVector& kv,
                                                   const size_t source_degree,
                                                   const size_t target_degree,
                                                   const double param_tol )
    {
        if( target_degree < source_degree ) throw std::invalid_argument( "Cannot degree elevate to a lower degree" );

        const size_t source_num_funcs = kv.size() - source_degree - 1;
        if( target_degree == source_degree )
        {
            Eigen::SparseMatrix<double> identity( source_num_funcs, source_num_funcs );
            identity.setIdentity();
            return identity;
        }

        if( kv.uniqueKnotMultiplicities().front().second != source_degree + 1 or
            kv.uniqueKnotMultiplicities().back().second != source_degree + 1 )
            throw std::runtime_error( "Degree elevation only supports open knot vectors" );

        const KnotVector target_kv = degreeElevatedKnotVector( kv, target_degree - source_degree );
        const size_t target_num_funcs = target_kv.size() - target_degree - 1;

        const std::vector<LocalExtraction> source_extractions = localExtractions( kv, source_degree );
        const std::vector<LocalExtraction> target_extractions = localExtractions( target_kv, target_degree );
        if( source_extractions.size() != target_extractions.size() )
            throw std::runtime_error( "Degree elevation changed the number of elements." );

        const Eigen::MatrixXd M = bernsteinDegreeElevationMatrix( source_degree, target_degree );
        const double duplicate_tol = duplicateTolerance( param_tol );

        std::map<std::pair<Eigen::Index, Eigen::Index>, std::pair<double, size_t>> entries;

        for( size_t elem_ii = 0; elem_ii < source_extractions.size(); elem_ii++ )
        {
            const LocalExtraction& source_local = source_extractions.at( elem_ii );
            const LocalExtraction& target_local = target_extractions.at( elem_ii );
            const Eigen::MatrixXd A = source_local.extraction * M;

            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver( target_local.extraction.transpose() );
            if( solver.rank() < target_local.extraction.rows() )
                throw std::runtime_error( "Singular local extraction operator in degree elevation." );

            const Eigen::MatrixXd local_op = solver.solve( A.transpose() ).transpose();

            for( Eigen::Index local_source = 0; local_source < local_op.rows(); local_source++ )
            {
                const Eigen::Index source_fid = source_local.connectivity.at( local_source ).id();
                for( Eigen::Index local_target = 0; local_target < local_op.cols(); local_target++ )
                {
                    const Eigen::Index target_fid = target_local.connectivity.at( local_target ).id();
                    const double value = local_op( local_source, local_target );
                    const auto key = std::make_pair( source_fid, target_fid );
                    auto [entry_it, inserted] = entries.emplace( key, std::make_pair( value, size_t{ 1 } ) );
                    if( not inserted )
                    {
                        const double average = entry_it->second.first / static_cast<double>( entry_it->second.second );
                        const double scale = std::max( { 1.0, std::abs( average ), std::abs( value ) } );
                        if( std::abs( average - value ) > duplicate_tol * scale )
                            throw std::runtime_error( "Inconsistent duplicate degree elevation transfer entry." );

                        entry_it->second.first += value;
                        entry_it->second.second++;
                    }
                }
            }
        }

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve( entries.size() );
        const double drop_tol = 100.0 * std::numeric_limits<double>::epsilon();
        for( const auto& [key, value_count] : entries )
        {
            const double value = value_count.first / static_cast<double>( value_count.second );
            if( std::abs( value ) > drop_tol )
                triplets.emplace_back( key.first, key.second, value );
        }

        Eigen::SparseMatrix<double> result( source_num_funcs, target_num_funcs );
        result.setFromTriplets( triplets.begin(), triplets.end() );
        return result;
    }
}
