#include <IndexOperations.hpp>
#include <numeric>

namespace util
{
    size_t flatten( const IndexVec& tp_indices, const IndexVec& lengths )
    {
        size_t multiplier = 1;
        size_t out = 0;
        for( size_t i = 0; i < tp_indices.size(); i++ )
        {
            if( tp_indices.at( i ) >= lengths.at( i ) )
                throw std::runtime_error( "TP index meets or exceeds given length" );

            out += multiplier * tp_indices.at( i );
            multiplier *= lengths.at( i );
        }
        return out;
    }

    IndexVec unflatten( const size_t index, const IndexVec& lengths )
    {
        IndexVec out;
        size_t divisor = 1;
        size_t mod_divisor = 1;
        for( size_t i = 0; i < lengths.size(); i++ )
        {
            mod_divisor *= lengths.at( i );
            out.push_back( ( index % mod_divisor ) / divisor );
            divisor *= lengths.at( i );
        }
        if( out.at( lengths.size() - 1 ) >= lengths.at( lengths.size() - 1 ) )
            throw std::runtime_error( "Flat index exceeds greatest index in TP space" );
        return out;
    }

    void iterateTensorProduct( const IndexVec& lengths, const std::function<void( const IndexVec& )>& callback )
    {
        const size_t num_ids = std::reduce( lengths.begin(), lengths.end(), 1, std::multiplies<>() );
        for( size_t id = 0; id < num_ids; id++ )
            callback( unflatten( id, lengths ) );
    }

    void iterateVTKTPOrdering( const IndexVec& lengths, const std::function<void( const size_t& )>& callback )
    {
        // There is probably a more unified way to write this...
        switch( lengths.size() )
        {
            case 1:
                callback( 0 );
                if( lengths.at( 0 ) > 1 ) callback( lengths.at( 0 ) - 1 );
                for( size_t i = 1; i < lengths.at( 0 ) - 1; i++ )
                    callback( i );
                break;
            case 2:
            {
                if( lengths.at( 0 ) < 2 or lengths.at( 1 ) < 2 )
                    throw std::runtime_error( "Cannot iterate with lengths < 2" ); //FIXME

                // Corners in CCW order
                callback( flatten( { 0, 0 }, lengths ) );
                callback( flatten( { lengths.at( 0 ) - 1, 0 }, lengths ) );
                callback( flatten( { lengths.at( 0 ) - 1, lengths.at( 1 ) - 1 }, lengths ) );
                callback( flatten( { 0, lengths.at( 1 ) - 1 }, lengths ) );

                // Edges in CCW order, iterated in positive s or t direction
                for( size_t i = 1; i < lengths.at( 0 ) - 1; i++ )
                    callback( flatten( {i, 0}, lengths ) );
                for( size_t i = 1; i < lengths.at( 1 ) - 1; i++ )
                    callback( flatten( {lengths.at( 0 ) - 1, i}, lengths ) );
                for( size_t i = 1; i < lengths.at( 0 ) - 1; i++ )
                    callback( flatten( {i, lengths.at( 1 ) - 1}, lengths ) );
                for( size_t i = 1; i < lengths.at( 1 ) - 1; i++ )
                    callback( flatten( {0, i}, lengths ) );

                // Interior in tensor product order
                iterateTensorProduct( { lengths.at( 0 ) - 2, lengths.at( 1 ) - 2 }, [&]( const IndexVec& interior_index ) {
                    callback( flatten( { interior_index.at( 0 ) + 1, interior_index.at( 1 ) + 1 }, lengths ) );
                } );
            }
            break;
            case 3:
            {
                if( lengths.at( 0 ) < 2 or lengths.at( 1 ) < 2 or lengths.at( 2 ) < 2 )
                    throw std::runtime_error( "Cannot iterate with lengths < 2" ); //FIXME

                // Corners in CCW order
                callback( flatten( { 0, 0, 0 }, lengths ) );
                callback( flatten( { lengths.at( 0 ) - 1, 0, 0 }, lengths ) );
                callback( flatten( { lengths.at( 0 ) - 1, lengths.at( 1 ) - 1, 0 }, lengths ) );
                callback( flatten( { 0, lengths.at( 1 ) - 1, 0 }, lengths ) );
                callback( flatten( { 0, 0, lengths.at( 2 ) - 1 }, lengths ) );
                callback( flatten( { lengths.at( 0 ) - 1, 0, lengths.at( 2 ) - 1 }, lengths ) );
                callback( flatten( { lengths.at( 0 ) - 1, lengths.at( 1 ) - 1, lengths.at( 2 ) - 1 }, lengths ) );
                callback( flatten( { 0, lengths.at( 1 ) - 1, lengths.at( 2 ) - 1 }, lengths ) );

                // Edges of u = 0 face in CCW order, iterated in positive s or t direction
                for( size_t i = 1; i < lengths.at( 0 ) - 1; i++ )
                    callback( flatten( {i, 0, 0}, lengths ) );
                for( size_t i = 1; i < lengths.at( 1 ) - 1; i++ )
                    callback( flatten( {lengths.at( 0 ) - 1, i, 0}, lengths ) );
                for( size_t i = 1; i < lengths.at( 0 ) - 1; i++ )
                    callback( flatten( {i, lengths.at( 1 ) - 1, 0}, lengths ) );
                for( size_t i = 1; i < lengths.at( 1 ) - 1; i++ )
                    callback( flatten( {0, i, 0}, lengths ) );
                // Edges of u = 1 face in CCW order, iterated in positive s or t direction
                for( size_t i = 1; i < lengths.at( 0 ) - 1; i++ )
                    callback( flatten( {i, 0, lengths.at( 2 ) - 1 }, lengths ) );
                for( size_t i = 1; i < lengths.at( 1 ) - 1; i++ )
                    callback( flatten( {lengths.at( 0 ) - 1, i, lengths.at( 2 ) - 1 }, lengths ) );
                for( size_t i = 1; i < lengths.at( 0 ) - 1; i++ )
                    callback( flatten( {i, lengths.at( 1 ) - 1, lengths.at( 2 ) - 1 }, lengths ) );
                for( size_t i = 1; i < lengths.at( 1 ) - 1; i++ )
                    callback( flatten( {0, i, lengths.at( 2 ) - 1 }, lengths ) );
                //Edges with changing u coordinate in TP order, iterated in positive u direction
                for( size_t i = 1; i < lengths.at( 2 ) - 1; i++ )
                    callback( flatten( { 0, 0, i }, lengths ) );
                for( size_t i = 1; i < lengths.at( 2 ) - 1; i++ )
                    callback( flatten( { lengths.at( 0 ) - 1, 0, i }, lengths ) );
                for( size_t i = 1; i < lengths.at( 2 ) - 1; i++ )
                    callback( flatten( { 0, lengths.at( 1 ) - 1, i }, lengths ) );
                for( size_t i = 1; i < lengths.at( 2 ) - 1; i++ )
                    callback( flatten( { lengths.at( 0 ) - 1, lengths.at( 1 ) - 1, i }, lengths ) );

                // s = 0 face
                iterateTensorProduct( { lengths.at( 1 ) - 2, lengths.at( 2 ) - 2 }, [&]( const IndexVec& interior_index ) {
                    callback( flatten( { 0, interior_index.at( 0 ) + 1, interior_index.at( 1 ) + 1 }, lengths ) );
                } );
                // s = 1 face
                iterateTensorProduct( { lengths.at( 1 ) - 2, lengths.at( 2 ) - 2 }, [&]( const IndexVec& interior_index ) {
                    callback( flatten( { lengths.at( 0 ) - 1, interior_index.at( 0 ) + 1, interior_index.at( 1 ) + 1 }, lengths ) );
                } );
                // t = 0 face
                iterateTensorProduct( { lengths.at( 0 ) - 2, lengths.at( 2 ) - 2 }, [&]( const IndexVec& interior_index ) {
                    callback( flatten( { interior_index.at( 0 ) + 1, 0, interior_index.at( 1 ) + 1 }, lengths ) );
                } );
                // t = 1 face
                iterateTensorProduct( { lengths.at( 0 ) - 2,  lengths.at( 2 ) - 2 }, [&]( const IndexVec& interior_index ) {
                    callback( flatten( { interior_index.at( 0 ) + 1, lengths.at( 1 ) - 1, interior_index.at( 1 ) + 1 }, lengths ) );
                } );
                // u = 0 face
                iterateTensorProduct( { lengths.at( 0 ) - 2, lengths.at( 1 ) - 2 }, [&]( const IndexVec& interior_index ) {
                    callback( flatten( { interior_index.at( 0 ) + 1, interior_index.at( 1 ) + 1, 0 }, lengths ) );
                } );
                // u = 1 face
                iterateTensorProduct( { lengths.at( 0 ) - 2, lengths.at( 1 ) - 2 }, [&]( const IndexVec& interior_index ) {
                    callback( flatten( { interior_index.at( 0 ) + 1, interior_index.at( 1 ) + 1, lengths.at( 2 ) - 1 }, lengths ) );
                } );

                // Interior in tensor product order
                iterateTensorProduct( { lengths.at( 0 ) - 2, lengths.at( 1 ) - 2, lengths.at( 2 ) - 2 }, [&]( const IndexVec& interior_index ) {
                    callback( flatten( { interior_index.at( 0 ) + 1, interior_index.at( 1 ) + 1, interior_index.at( 2 ) + 1 }, lengths ) );
                } );
            }
            break;
            default:
            break;
        }
    }

    Eigen::MatrixXd tensorProduct( const SmallVector<Eigen::VectorXd, 3>& vecs )
    {
        const SmallVector<size_t, 3> lengths = std::accumulate( vecs.begin(), vecs.end(), SmallVector<size_t, 3>(), [&]( SmallVector<size_t, 3> accum, const Eigen::VectorXd& vec ) {
            accum.push_back( vec.size() );
            return accum;
        } );
        const size_t rows = std::accumulate( lengths.begin(), lengths.end(), 1, []( const size_t a, const size_t b ) { return a * b; } );

        Eigen::MatrixXd result( rows, vecs.size() );
        iterateTensorProduct( lengths, [&]( const IndexVec& ids ){
            const size_t row = flatten( ids, lengths );
            for( size_t i = 0; i < vecs.size(); i++ )
                result( row, i ) = vecs.at( i )( ids.at( i ) );
        } );

        return result;
    }
}