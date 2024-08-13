#include <TPCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>
#include <utility>

using namespace topology;

// The permutation of TPDartPos to perform a phi(1) or phi(-1) in 2d
inline TPCombinatorialMap::TPDartPos phi1s_2d( const TPCombinatorialMap::TPDartPos& pos, const int i )
{
    return TPCombinatorialMap::TPDartPos( ( std::to_underlying( pos ) + i ) % 4 );
}

// The permutation of TPDartPos to perform a phi(2) in 2d
inline TPCombinatorialMap::TPDartPos phi2s_2d( const TPCombinatorialMap::TPDartPos& pos )
{
    switch( pos )
    {
        case TPCombinatorialMap::TPDartPos::DartPos0: return TPCombinatorialMap::TPDartPos::DartPos2;
        case TPCombinatorialMap::TPDartPos::DartPos1: return TPCombinatorialMap::TPDartPos::DartPos3;
        case TPCombinatorialMap::TPDartPos::DartPos2: return TPCombinatorialMap::TPDartPos::DartPos0;
        case TPCombinatorialMap::TPDartPos::DartPos3: return TPCombinatorialMap::TPDartPos::DartPos1;
        default: throw std::runtime_error( "Cannot have DartPos > 3 in 2d" );
    }
}

TPCombinatorialMap::TPCombinatorialMap( const std::shared_ptr<const CombinatorialMap>& source,
                                        const std::shared_ptr<const CombinatorialMap1d>& line )
    : mSource( source ), mLine( line )
{}

std::optional<Dart> TPCombinatorialMap::phi( const int i, const Dart& d ) const
{
    const auto [source_dart, line_dart, pos] = unflatten( d );
    if( dim() == 2 )
    {
        if( i == 1 or i == -1 ) return flatten( source_dart, line_dart, phi1s_2d( pos, i ) );
        if( i == 2 )
        {
            switch( pos )
            {
                case TPDartPos::DartPos0:
                    return topology::phi( *mLine, -1, line_dart ).transform( [&]( const Dart& new_line_dart ){
                        return flatten( source_dart, new_line_dart, phi2s_2d( pos ) );
                    } );
                case TPDartPos::DartPos1:
                    return topology::phi( *mSource, 1, source_dart ).transform( [&]( const Dart& new_source_dart ){
                        return flatten( new_source_dart, line_dart, phi2s_2d( pos ) );
                    } );
                case TPDartPos::DartPos2:
                    return topology::phi( *mLine, 1, line_dart ).transform( [&]( const Dart& new_line_dart ){
                        return flatten( source_dart, new_line_dart, phi2s_2d( pos ) );
                    } );
                case TPDartPos::DartPos3:
                    return topology::phi( *mSource, -1, source_dart ).transform( [&]( const Dart& new_source_dart ){
                        return flatten( new_source_dart, line_dart, phi2s_2d( pos ) );
                    } );
                default: throw std::runtime_error( "Cannot have DartPos > 3 in 2d" );
            }
        }
    }
    else if( dim() == 3 )
    {
        switch( pos )
        {
            case TPDartPos::DartPos0:
            {
                switch( i )
                {
                    case -1:
                    case 1:
                        return topology::phi( *mSource, i, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flatten( new_source_dart, line_dart, pos );
                        } );
                    case 2:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos1 );
                    case 3:
                        return topology::phi( *mLine, -1, line_dart ).transform( [&]( const Dart& new_line_dart ){
                            return flatten( source_dart, new_line_dart, TPDartPos::DartPos5 );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos1:
            {
                switch( i )
                {
                    case -1:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos4 );
                    case 1:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos2 );
                    case 2:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos0 );
                    case 3:
                        return topology::phi( *mSource, 2, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flatten( new_source_dart, line_dart, pos );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos2:
            {
                switch( i )
                {
                    case -1:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos1 );
                    case 1:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos3 );
                    case 2:
                        return topology::phi( *mSource, -1, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flatten( new_source_dart, line_dart, TPDartPos::DartPos4 );
                        } );
                    case 3:
                        return topology::phi( *mSource, 2, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flatten( new_source_dart, line_dart, TPDartPos::DartPos4 );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos3:
            {
                switch( i )
                {
                    case -1:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos2 );
                    case 1:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos4 );
                    case 2:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos5 );
                    case 3:
                        return topology::phi( *mSource, 2, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flatten( new_source_dart, line_dart, pos );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos4:
            {
                switch( i )
                {
                    case -1:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos3 );
                    case 1:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos1 );
                    case 2:
                        return topology::phi( *mSource, 1, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flatten( new_source_dart, line_dart, TPDartPos::DartPos2 );
                        } );
                    case 3:
                        return topology::phi( *mSource, 2, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flatten( new_source_dart, line_dart, TPDartPos::DartPos2 );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos5:
            {
                switch( i )
                {
                    case -1:
                    case 1:
                        return topology::phi( *mSource, -1 * i, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flatten( new_source_dart, line_dart, pos );
                        } );
                    case 2:
                        return flatten( source_dart, line_dart, TPDartPos::DartPos3 );
                    case 3:
                        return topology::phi( *mLine, 1, line_dart ).transform( [&]( const Dart& new_line_dart ){
                            return flatten( source_dart, new_line_dart, TPDartPos::DartPos0 );
                        } );
                    default: return std::nullopt;
                }
            }
        }
    }

    return std::nullopt;
}

Dart::IndexType TPCombinatorialMap::maxDartId() const
{
    return ( mSource->maxDartId() + 1 ) * ( mLine->maxDartId() + 1 ) * 2 * dim() - 1;
}

uint TPCombinatorialMap::dim() const
{
    return mSource->dim() + 1;
}

bool TPCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    return topology::iterateDartsWhile( *mSource, [&]( const Dart& source_dart ) {
        return topology::iterateDartsWhile( *mLine, [&]( const Dart& line_dart ) {
            for( Dart::IndexType pos = 0; pos < dartsPerSourceDart(); pos++ )
            {
                if( not callback( flatten( source_dart, line_dart, TPDartPos( pos ) ) ) ) return false;
            }
            return true;
        } );
    } );
}

bool TPCombinatorialMap::iterateCellsWhile( const uint cell_dim,
                                const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim > dim() ) return true;
    if( cell_dim == dim() )
    {
        return topology::iterateDartsWhile( *mSource, [&]( const Dart& source_dart ) {
            return topology::iterateDartsWhile( *mLine, [&]( const Dart& line_dart ) {
                return callback( Cell( flatten( source_dart, line_dart, TPDartPos::DartPos0 ), cell_dim ) );
            } );
        } );
    }
    else
    {
        GlobalCellMarker m( *this, cell_dim );
        return iterateDartsWhile( [&]( const Dart& d ){
            const Cell c( d, cell_dim );
            if( not m.isMarked( c ) )
            {
                m.mark( *this, c );
                if( not callback( c ) ) return false;
            }
            return true;
        } );
    }
    return true;
}

std::optional<IndexingFunc> TPCombinatorialMap::indexing( const uint ) const
{
    return std::nullopt;
}

std::optional<size_t> TPCombinatorialMap::cellCount( const uint cell_dim ) const
{
    if( cell_dim == dim() )
    {
        return mSource->cellCount( dim() - 1 ).and_then( [&]( const size_t a ) {
            return mLine->cellCount( 1 ).transform( [&a]( const size_t b ) { return a * b; } );
        } );
    }

    // TODO: There is a way to calculate this from the underlying cell counts.

    return std::nullopt;
}

Dart TPCombinatorialMap::flatten( const Dart& source_dart, const Dart& line_dart, const TPDartPos& pos ) const
{
    return Dart( std::to_underlying( pos ) + source_dart.id() * dartsPerSourceDart() +
                 line_dart.id() * ( mSource->maxDartId() + 1 ) * dartsPerSourceDart() );
}

std::tuple<Dart, Dart, TPCombinatorialMap::TPDartPos> TPCombinatorialMap::unflatten( const Dart& d ) const
{
    return { Dart( ( d.id() / dartsPerSourceDart() ) % ( mSource->maxDartId() + 1 ) ),
             Dart( d.id() / ( ( mSource->maxDartId() + 1 ) * dartsPerSourceDart() ) ),
             TPDartPos( d.id() % dartsPerSourceDart() ) };
}

namespace topology
{
    std::vector<std::shared_ptr<const CombinatorialMap1d>> tensorProductComponentCMaps( const TPCombinatorialMap& tp_map )
    {
        const size_t dim = tp_map.dim();
        std::vector<std::shared_ptr<const CombinatorialMap1d>> out;
        out.reserve( dim );

        if( dim == 3 )
        {
            const std::shared_ptr<const TPCombinatorialMap> source_primal =
                std::dynamic_pointer_cast<const TPCombinatorialMap>( tp_map.sourceCMapPtr() );
            if( source_primal.get() == nullptr ) return {};
            out.push_back( std::dynamic_pointer_cast<const CombinatorialMap1d>( source_primal->sourceCMapPtr() ) );
            if( out.back().get() == nullptr ) return {};
            out.push_back( source_primal->lineCMapPtr() );
        }
        else
        {
            out.push_back( std::dynamic_pointer_cast<const CombinatorialMap1d>( tp_map.sourceCMapPtr() ) );
            if( out.back().get() == nullptr ) return {};
        }
        out.push_back( tp_map.lineCMapPtr() );
        return out;
    }
}