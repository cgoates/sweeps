#pragma once
#include <vector>
#include <stdexcept>
#include <optional>

namespace util
{
    class UnionFind
    {
        public:
        UnionFind( const size_t n ) : parent( n ), size( n, 1 ), is_opposite_parent( n, false ), num_sets( n )
        {
            for( size_t i = 0; i < n; i++ )
            {
                parent.at( i ) = i;
            }
        }

        void unite( const size_t a, const size_t b, const bool aligned )
        {
            const size_t root_a = find( a );
            const size_t root_b = find( b );

            const bool a_is_opposite = is_opposite_parent.at( a );
            const bool b_is_opposite = is_opposite_parent.at( b );

            if( root_a == root_b )
            {
                // Already in the same set, check for consistency
                if( ( a_is_opposite == b_is_opposite ) != aligned )
                {
                    throw std::runtime_error( "UnionFind: trying to unite two elements with inconsistent orientation!" );
                }
                return;
            }

            if( --num_sets == 0 )
            {
                throw std::runtime_error( "UnionFind: trying to unite more sets than exist!" );
            }

            const bool parity_between_roots = ( a_is_opposite != b_is_opposite ) == aligned;

            //Union by size
            if( size.at( root_a ) < size.at( root_b ) )
            {
                parent.at( root_a ) = root_b;
                is_opposite_parent.at( root_a ) = parity_between_roots;
                size.at( root_b ) += size.at( root_a );
            }
            else
            {
                parent.at( root_b ) = root_a;
                is_opposite_parent.at( root_b ) = parity_between_roots;
                size.at( root_a ) += size.at( root_b );
            }
        }

        std::pair<size_t, bool> findWithOrientation( const size_t a )
        {
            const size_t root = find( a );
            return { root, is_opposite_parent.at( a ) };
        }

        size_t numSets() const
        {
            return num_sets;
        }

        private:

        size_t find( const size_t a )
        {
            if( parent.at( a ) == a ) return a;

            const size_t p = parent.at( a );
            const size_t root = find( p );
            is_opposite_parent.at( a ) = ( is_opposite_parent.at( a ) != is_opposite_parent.at( p ) );
            parent.at( a ) = root;
            return root;
        }

        std::vector<size_t> parent;
        std::vector<size_t> size;
        std::vector<bool> is_opposite_parent;
        size_t num_sets;
    };
}