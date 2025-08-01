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

Dart flattenInternal( const Dart& source_dart, const Dart& line_dart, const TPCombinatorialMap::TPDartPos& pos, size_t darts_per_source_dart, Dart::IndexType source_max_dart )
{
    return Dart( std::to_underlying( pos ) + source_dart.id() * darts_per_source_dart +
                 line_dart.id() * ( source_max_dart + 1 ) * darts_per_source_dart );
}

std::tuple<Dart, Dart, TPCombinatorialMap::TPDartPos> unflattenInternal( const Dart& d, size_t darts_per_source_dart, Dart::IndexType source_max_dart )
{
    return { Dart( ( d.id() / darts_per_source_dart ) % ( source_max_dart + 1 ) ),
             Dart( d.id() / ( ( source_max_dart + 1 ) * darts_per_source_dart ) ),
             TPCombinatorialMap::TPDartPos( d.id() % darts_per_source_dart ) };
}

std::optional<Dart> TPCombinatorialMap::phi( const int i, const Dart& d ) const
{
    const size_t darts_per_source_dart = dartsPerSourceDart();
    const Dart::IndexType source_max_dart = mSource->maxDartId();
    const uint dim = this->dim();
    const auto [source_dart, line_dart, pos] = unflattenInternal( d, darts_per_source_dart, source_max_dart );
    if( dim == 2 )
    {
        if( i == 1 or i == -1 ) return flattenInternal( source_dart, line_dart, phi1s_2d( pos, i ), darts_per_source_dart, source_max_dart );
        if( i == 2 )
        {
            switch( pos )
            {
                case TPDartPos::DartPos0:
                    return topology::phi( *mLine, -1, line_dart ).transform( [&]( const Dart& new_line_dart ){
                        return flattenInternal( source_dart, new_line_dart, phi2s_2d( pos ), darts_per_source_dart, source_max_dart );
                    } );
                case TPDartPos::DartPos1:
                    return topology::phi( *mSource, 1, source_dart ).transform( [&]( const Dart& new_source_dart ){
                        return flattenInternal( new_source_dart, line_dart, phi2s_2d( pos ), darts_per_source_dart, source_max_dart );
                    } );
                case TPDartPos::DartPos2:
                    return topology::phi( *mLine, 1, line_dart ).transform( [&]( const Dart& new_line_dart ){
                        return flattenInternal( source_dart, new_line_dart, phi2s_2d( pos ), darts_per_source_dart, source_max_dart );
                    } );
                case TPDartPos::DartPos3:
                    return topology::phi( *mSource, -1, source_dart ).transform( [&]( const Dart& new_source_dart ){
                        return flattenInternal( new_source_dart, line_dart, phi2s_2d( pos ), darts_per_source_dart, source_max_dart );
                    } );
                default: throw std::runtime_error( "Cannot have DartPos > 3 in 2d" );
            }
        }
    }
    else if( dim == 3 )
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
                            return flattenInternal( new_source_dart, line_dart, pos, darts_per_source_dart, source_max_dart );
                        } );
                    case 2:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos1, darts_per_source_dart, source_max_dart );
                    case 3:
                        return topology::phi( *mLine, -1, line_dart ).transform( [&]( const Dart& new_line_dart ){
                            return flattenInternal( source_dart, new_line_dart, TPDartPos::DartPos5, darts_per_source_dart, source_max_dart );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos1:
            {
                switch( i )
                {
                    case -1:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos4, darts_per_source_dart, source_max_dart );
                    case 1:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos2, darts_per_source_dart, source_max_dart );
                    case 2:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos0, darts_per_source_dart, source_max_dart );
                    case 3:
                        return topology::phi( *mSource, 2, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flattenInternal( new_source_dart, line_dart, pos, darts_per_source_dart, source_max_dart );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos2:
            {
                switch( i )
                {
                    case -1:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos1, darts_per_source_dart, source_max_dart );
                    case 1:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos3, darts_per_source_dart, source_max_dart );
                    case 2:
                        return topology::phi( *mSource, -1, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flattenInternal( new_source_dart, line_dart, TPDartPos::DartPos4, darts_per_source_dart, source_max_dart );
                        } );
                    case 3:
                        return topology::phi( *mSource, 2, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flattenInternal( new_source_dart, line_dart, TPDartPos::DartPos4, darts_per_source_dart, source_max_dart );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos3:
            {
                switch( i )
                {
                    case -1:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos2, darts_per_source_dart, source_max_dart );
                    case 1:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos4, darts_per_source_dart, source_max_dart );
                    case 2:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos5, darts_per_source_dart, source_max_dart );
                    case 3:
                        return topology::phi( *mSource, 2, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flattenInternal( new_source_dart, line_dart, pos, darts_per_source_dart, source_max_dart );
                        } );
                    default: return std::nullopt;
                }
            }
            case TPDartPos::DartPos4:
            {
                switch( i )
                {
                    case -1:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos3, darts_per_source_dart, source_max_dart );
                    case 1:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos1, darts_per_source_dart, source_max_dart );
                    case 2:
                        return topology::phi( *mSource, 1, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flattenInternal( new_source_dart, line_dart, TPDartPos::DartPos2, darts_per_source_dart, source_max_dart );
                        } );
                    case 3:
                        return topology::phi( *mSource, 2, source_dart ).transform( [&]( const Dart& new_source_dart ){
                            return flattenInternal( new_source_dart, line_dart, TPDartPos::DartPos2, darts_per_source_dart, source_max_dart );
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
                            return flattenInternal( new_source_dart, line_dart, pos, darts_per_source_dart, source_max_dart );
                        } );
                    case 2:
                        return flattenInternal( source_dart, line_dart, TPDartPos::DartPos3, darts_per_source_dart, source_max_dart );
                    case 3:
                        return topology::phi( *mLine, 1, line_dart ).transform( [&]( const Dart& new_line_dart ){
                            return flattenInternal( source_dart, new_line_dart, TPDartPos::DartPos0, darts_per_source_dart, source_max_dart );
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
    const size_t darts_per_source_dart = dartsPerSourceDart();
    const Dart::IndexType source_max_dart = mSource->maxDartId();
    return topology::iterateDartsWhile( *mSource, [&]( const Dart& source_dart ) {
        return topology::iterateDartsWhile( *mLine, [&]( const Dart& line_dart ) {
            for( Dart::IndexType pos = 0; pos < dartsPerSourceDart(); pos++ )
            {
                if( not callback( flattenInternal( source_dart, line_dart, TPDartPos( pos ), darts_per_source_dart, source_max_dart ) ) ) return false;
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
        const size_t darts_per_source_dart = dartsPerSourceDart();
        const Dart::IndexType source_max_dart = mSource->maxDartId();
        return topology::iterateCellsWhile( *mSource, cell_dim - 1, [&]( const Cell& source_cell ) {
            return topology::iterateDartsWhile( *mLine, [&]( const Dart& line_dart ) {
                return callback( Cell( flattenInternal( source_cell.dart(), line_dart, TPDartPos::DartPos0, darts_per_source_dart, source_max_dart ), cell_dim ) );
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
    return flattenInternal( source_dart, line_dart, pos, dartsPerSourceDart(), mSource->maxDartId() );
}

std::tuple<Dart, Dart, TPCombinatorialMap::TPDartPos> TPCombinatorialMap::unflatten( const Dart& d ) const
{
    return unflattenInternal( d, dartsPerSourceDart(), mSource->maxDartId() );
}

namespace topology
{
    SmallVector<std::shared_ptr<const CombinatorialMap1d>, 3> tensorProductComponentCMaps( const TPCombinatorialMap& tp_map )
    {
        const size_t dim = tp_map.dim();
        SmallVector<std::shared_ptr<const CombinatorialMap1d>, 3> out;

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

    FullyUnflattenedDart unflattenFull( const TPCombinatorialMap& cmap, const Dart& d )
    {
        const size_t dim = cmap.dim();
        FullyUnflattenedDart out;

        const auto unflat = cmap.unflatten( d );

        if( dim == 3 )
        {
            const std::shared_ptr<const TPCombinatorialMap> source =
                std::dynamic_pointer_cast<const TPCombinatorialMap>( cmap.sourceCMapPtr() );
            if( source.get() == nullptr )
                out.unflat_darts.push_back( std::get<0>( unflat ) );
            else
            {
                const auto unflat2 = source->unflatten( std::get<0>( unflat ) );
                out.unflat_darts.push_back( std::get<0>( unflat2 ) );
                out.unflat_darts.push_back( std::get<1>( unflat2 ) );
                out.dart_pos.push_back( std::get<2>( unflat2 ) );
            }
        }
        else
        {
            out.unflat_darts.push_back( std::get<0>( unflat ) );
        }
        out.unflat_darts.push_back( std::get<1>( unflat ) );
        out.dart_pos.push_back( std::get<2>( unflat ) );
        return out;
    }

    TPCombinatorialMap tensorProductCMapFromComponents( const SmallVector<std::shared_ptr<const CombinatorialMap1d>, 3>& components )
    {
        if( components.size() == 2 )
        {
            return TPCombinatorialMap( components.at( 0 ), components.at( 1 ) );
        }
        else if( components.size() == 3 )
        {
            const auto source = std::make_shared<const TPCombinatorialMap>( components.at( 0 ), components.at( 1 ) );
            return TPCombinatorialMap( source, components.at( 2 ) );
        }
        else
        {
            throw std::invalid_argument( "TPCombinatorialMap must have 2 or 3 components" );
        }
    }

    Dart flattenFull( const TPCombinatorialMap& cmap, const FullyUnflattenedDart& unflat_darts )
    {
        if( cmap.dim() == 3 )
        {
            const std::shared_ptr<const TPCombinatorialMap> source =
                std::dynamic_pointer_cast<const TPCombinatorialMap>( cmap.sourceCMapPtr() );
            if( source.get() == nullptr )
                throw std::runtime_error( "Cannot flattenFull on 3d TP cmap with non TP source." );

            const Dart source_d =
                source->flatten( unflat_darts.unflat_darts.at( 0 ), unflat_darts.unflat_darts.at( 1 ), unflat_darts.dart_pos.at( 0 ) );
            return cmap.flatten( source_d, unflat_darts.unflat_darts.at( 2 ), unflat_darts.dart_pos.at( 1 ) );
        }
        else
        {
            return cmap.flatten( unflat_darts.unflat_darts.at( 0 ), unflat_darts.unflat_darts.at( 1 ), unflat_darts.dart_pos.at( 0 ) );
        }
    }

    std::pair<CellOrEndVertex, CellOrEndVertex> unflattenCell( const TPCombinatorialMap& cmap, const Cell& c )
    {
        const size_t cmap_dim = cmap.dim();
        if( c.dim() > cmap_dim )
            throw std::invalid_argument( "Cannot unflatten cell with dimension greater than the cmap" );

        const auto darts = cmap.unflatten( c.dart() );
        const auto& tp_pos = std::get<2>( darts );
        const auto tp_pos_int = std::to_underlying( tp_pos );
        const auto& line_d = std::get<1>( darts );
        const auto& source_d = std::get<0>( darts );

        if( c.dim() == cmap_dim )
        {
            return { Cell( source_d, cmap_dim - 1 ), Edge( line_d ) };
        }
        
        const auto opp_vert = []( const CombinatorialMap& cmap, const Dart& d ) -> CellOrEndVertex {
            return phi( cmap, 1, d ).and_then( []( const Dart& a ) -> std::optional<Vertex> {
                return Vertex( a );
            } );
        };
        // If it's the line and you know it should be a vert, just call this.
        const auto line_vert = [&]() -> CellOrEndVertex {
            if( tp_pos_int >= cmap_dim )
                return opp_vert( cmap.lineCMap(), line_d );
            else
                return Vertex( line_d );
        };

        // If it's the source and you know it should be a vert, just call this.
        const auto source_vert = [&]() -> CellOrEndVertex {
            if( cmap_dim == 3 )
            {
                switch( tp_pos )
                {
                    case topology::TPCombinatorialMap::TPDartPos::DartPos0:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos2:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos3:
                        return Vertex( source_d );
                    default: return opp_vert( cmap.sourceCMap(), source_d );
                }
            }
            else
            {
                switch( tp_pos )
                {
                    case topology::TPCombinatorialMap::TPDartPos::DartPos0:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos3:
                        return Vertex( source_d );
                    default:
                        return opp_vert( cmap.sourceCMap(), source_d );
                }
            }
        };

        const auto line_edge = [&line_d](){
            return Edge( line_d );
        };

        const auto source_edge = [&source_d](){
            return Edge( source_d );
        };

        const auto source_face = [&source_d](){
            return Face( source_d );
        };

        if( c.dim() == 0 )
        {
            return { source_vert(), line_vert() };
        }
        else if( cmap_dim == 2 )
        {
            if( tp_pos_int % 2 == 0 )
                return { source_edge(), line_vert() };
            else
                return { source_vert(), line_edge() };
        }
        else//cmap.dim() == 3
        {
            if( c.dim() == 1 )
            {
                switch( tp_pos )
                {
                    case topology::TPCombinatorialMap::TPDartPos::DartPos2:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos4:
                        return { source_vert(), line_edge() };
                    default:
                        return { source_edge(), line_vert() };
                }
            }
            else
            {
                switch( tp_pos )
                {
                    case topology::TPCombinatorialMap::TPDartPos::DartPos0:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos5:
                        return { source_face(), line_vert() };
                    default:
                        return { source_edge(), line_edge() };
                }
            }
        }
    }
}