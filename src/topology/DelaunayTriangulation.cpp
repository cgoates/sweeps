#include <DelaunayTriangulation.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>

using namespace topology;

double ncos( const double adj0, const double adj1, const double opp )
{
    return -( adj0 * adj0 + adj1 * adj1 - opp * opp ) / ( adj0 * adj0 );
}

bool shouldConnectFirstThree( const std::array<Eigen::Vector3d, 4>& points )
{
    const double a = ( points[1] - points[0] ).norm();
    const double b = ( points[2] - points[1] ).norm();
    const double c = ( points[3] - points[2] ).norm();
    const double d = ( points[0] - points[3] ).norm();
    const double e = ( points[2] - points[0] ).norm();
    const double f = ( points[3] - points[1] ).norm();

    return std::min( { ncos( e, a, b ), ncos( e, b, a ), ncos( e, c, d ), ncos( e, d, c ) } ) >
            std::min( { ncos( f, a, d ), ncos( f, d, a ), ncos( f, c, b ), ncos( f, b, c ) } );
}


DelaunayTriangulation::DelaunayTriangulation( const CombinatorialMap& base, const VertexPositionsFunc& vert_positions )
    : mBaseMap( base ), mLowerBound( base.maxDartId() + base.maxDartId() % 2 ), mMaxDartId( mLowerBound )
{
    //TODO: Add the extra edges
    const auto split_face = [&]( const Face& f ) {
        const Dart new_d1( ++mMaxDartId );
        const Dart new_d2( ++mMaxDartId );

        const Dart old_phi_1 = topology::phi( base, -1, f.dart() ).value();
        const Dart old_phi1 = topology::phi( base, 1, f.dart() ).value();
        const Dart old_phi11 = topology::phi( base, {1, 1}, f.dart() ).value();

        mAlteredPhi1s.emplace( new_d1, f.dart() );
        mAlteredPhi_1s.emplace( f.dart(), new_d1 );

        mAlteredPhi1s.emplace( old_phi1, new_d1 );
        mAlteredPhi_1s.emplace( new_d1, old_phi1 );

        mAlteredPhi1s.emplace( old_phi_1, new_d2 );
        mAlteredPhi_1s.emplace( new_d2, old_phi_1 );

        mAlteredPhi1s.emplace( new_d2, old_phi11 );
        mAlteredPhi_1s.emplace( old_phi11, new_d2 );
    };

    topology::iterateCellsWhile( base, 2, [&]( const Face& f ) {
        size_t n_verts = 0;
        std::array<Eigen::Vector3d, 4> points;
        Dart d = f.dart();
        do
        {
            points[ n_verts ] = vert_positions( Vertex( d ) );
            n_verts++;
            d = topology::phi( base, 1, d ).value();
        } while( d != f.dart() );

        if( n_verts == 3 ) return true;
        if( n_verts > 4 ) throw std::runtime_error( "DelaunayTriangulation only handles quads" );
        if( shouldConnectFirstThree( points ) )
        {
            split_face( f );
        }
        else
        {
            split_face( topology::phi( base, 1, f.dart() ).value() );
        }
        return true;
    } );
}

std::optional<Dart> DelaunayTriangulation::phi( const int i, const Dart& d ) const
{
    if( i == 1 )
    {
        const auto it = mAlteredPhi1s.find( d );
        if( it == mAlteredPhi1s.end() )
        {
            if( d.id() > mLowerBound ) throw std::runtime_error( "Missing phi1 of dart " + std::to_string( d.id() ) );
            else return topology::phi( mBaseMap, i, d );
        }
        else
            return it->second;
    }
    else if( i == -1 )
    {
        const auto it = mAlteredPhi_1s.find( d );
        if( it == mAlteredPhi_1s.end() )
        {
            if( d.id() > mLowerBound ) throw std::runtime_error( "Missing phi-1 of dart " + std::to_string( d.id() ) );
            else return topology::phi( mBaseMap, i, d );
        }
        else
            return it->second;
    }
    else if( i == 2 )
    {
        if( d.id() > mLowerBound )
            return Dart( d.id() + ( d.id() % 2 == 0 ? -1 : 1 ) );
        else
            return topology::phi( mBaseMap, i, d );
    }
    else
        throw std::runtime_error( "Bad phi value" );
}

bool DelaunayTriangulation::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    const bool keep_iterating = topology::iterateDartsWhile( mBaseMap, callback );
    if( keep_iterating )
    {
        for( Dart::IndexType i = mLowerBound + 1; i <= mMaxDartId; i++ )
        {
            if( not callback( Dart( i ) ) ) return false;
        }
    }
    return keep_iterating;
}

bool DelaunayTriangulation::iterateCellsWhile( const uint cell_dim,
                                               const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim == 0 )
    {
        return topology::iterateCellsWhile( mBaseMap, cell_dim, callback );
    }
    else if( cell_dim == 1 )
    {
        const bool keep_iterating = topology::iterateCellsWhile( mBaseMap, cell_dim, callback );
        if( keep_iterating )
        {
            for( Dart::IndexType i = mLowerBound + 1; i <= mMaxDartId; i = i + 2 )
            {
                if( not callback( Edge( Dart( i ) ) ) ) return false;
            }
        }
        return keep_iterating;
    }
    else
    {
        GlobalCellMarker m( *this, cell_dim );
        const bool keep_iterating = topology::iterateCellsWhile( mBaseMap, cell_dim, [&]( const Face& f ) {
            m.mark( *this, f );
            return callback( f );
        } );
        if( keep_iterating )
        {
            for( Dart::IndexType i = mLowerBound + 1; i <= mMaxDartId; i++ )
            {
                const Face f( Dart{ i } );
                if( not m.isMarked( f ) )
                {
                    m.mark( *this, f );
                    if( not callback( f ) ) return false;
                }
            }
        }
        return keep_iterating;
    }
}

std::optional<IndexingFunc> DelaunayTriangulation::indexing( const uint cell_dim ) const
{
    if( cell_dim == 0 )
    {
        return mBaseMap.indexing( cell_dim )
            .and_then( [this]( const IndexingFunc& underlying_indexing ) -> std::optional<IndexingFunc> {
                return [this, underlying_indexing]( const Vertex& v ) {
                    std::optional<size_t> out = std::nullopt;

                    // Find a dart on the vertex from the underlying cmap
                    iterateDartsOfCell( *this, v, [&]( const Dart& d ) {
                        if( d.id() <= mLowerBound )
                        {
                            // and use the underlying index
                            out.emplace( underlying_indexing( Vertex( d ) ) );
                            return false;
                        }
                        return true;
                    } );
                    if( not out.has_value() )
                        throw std::runtime_error( "Vertex has no darts from underlying map!" );
                    return out.value();
                };
            } );
    }
    // NOTE: Could add edges if we store the max underlying edge id
    return std::nullopt;
}

std::optional<size_t> DelaunayTriangulation::cellCount( const uint cell_dim ) const
{
    if( cell_dim == 0 ) return mBaseMap.cellCount( cell_dim );
    else return topology::cellCount( mBaseMap, cell_dim ) + ( mMaxDartId - mLowerBound ) / 2;
}

namespace topology
{
    VertexPositionsFunc delaunayTriangulationVertexPositions(
        const DelaunayTriangulation& dtri, const VertexPositionsFunc& underlying_positions )
    {
        Dart::IndexType largest_underlying_dart_id = dtri.baseMap().maxDartId();
        return [&dtri, largest_underlying_dart_id, underlying_positions]( const Vertex& v ) -> Eigen::Vector3d {
            Eigen::Vector3d out;
            const bool found_it = not iterateDartsOfCell( dtri, v, [&]( const Dart& d ) {
                if( d.id() <= largest_underlying_dart_id )
                {
                    out = underlying_positions( Vertex( d ) );
                    return false;
                }
                return true;
            } );
            if( not found_it ) throw std::runtime_error( "No dart in vertex from underlying cmap!" );
            return out;
        };
    }
}