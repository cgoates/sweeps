#include <CombinatorialMap1d.hpp>

using namespace topology;

CombinatorialMap1d::CombinatorialMap1d( const Dart::IndexType n_cells, const bool periodic )
    : mNumCells( n_cells ), mPeriodic( periodic )
{}

std::optional<Dart> CombinatorialMap1d::phi( const int i, const Dart& d ) const
{
    if( i > 1 ) return std::nullopt;
    if( i == 1 and d == maxDartId() )
    {
        if( mPeriodic )
            return Dart( 0 );
        else
            return std::nullopt;
    }
    if( i == -1 and d == 0 )
    {
        if( mPeriodic )
            return Dart( maxDartId() );
        else
            return std::nullopt;
    }

    return Dart( d.id() + i );
}

Dart::IndexType CombinatorialMap1d::maxDartId() const
{
    return mNumCells - 1;
}

uint CombinatorialMap1d::dim() const
{
    return 1;
}

bool CombinatorialMap1d::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    for( Dart::IndexType d = 0; d < mNumCells; d++ )
    {
        if( not callback( Dart( d ) ) ) return false;
    }
    return true;
}

bool CombinatorialMap1d::iterateCellsWhile( const uint cell_dim,
                                            const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim > 1 ) return true;

    return iterateDartsWhile( [&]( const Dart& d ) { return callback( Cell( d, cell_dim ) ); } );
}

std::optional<IndexingFunc> CombinatorialMap1d::indexing( const uint cell_dim ) const
{
    if( cell_dim > 1 ) return std::nullopt;
    return []( const Cell& c ) { return c.dart().id(); };
}

std::optional<size_t> CombinatorialMap1d::cellCount( const uint cell_dim ) const
{
    if( cell_dim > 1 )
        return 0;
    else
        return mNumCells;
}