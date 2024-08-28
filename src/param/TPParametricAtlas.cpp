#include <TPParametricAtlas.hpp>
#include <utility>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>

namespace param
{
    TPParametricAtlas::TPParametricAtlas( const std::shared_ptr<const topology::TPCombinatorialMap>& cmap,
                                          const std::shared_ptr<const ParametricAtlas>& source_atlas,
                                          const std::shared_ptr<const ParametricAtlas1d>& line_atlas )
        : mMap( cmap ), mSourceParam( source_atlas ), mLineParam( line_atlas )
    {
        if( &cmap->sourceCMap() != &source_atlas->cmap() or &cmap->lineCMap() != &line_atlas->cmap() )
            throw std::runtime_error( "Incompatible cmap and atlases provided to TPParametricAtlas" );
    }

    const ParentDomain TPParametricAtlas::parentDomain( const topology::Cell& c ) const
    {
        const auto [source_dart, line_dart, tp_pos] = mMap->unflatten( c.dart() );
        const ParentDomain source_domain =
            mSourceParam->parentDomain( topology::Cell( source_dart, mSourceParam->cmap().dim() ) );
        const ParentDomain line_domain =
            mLineParam->parentDomain( topology::Cell( line_dart, mLineParam->cmap().dim() ) );
        return tensorProduct( source_domain, line_domain );
    }

    ParentPoint TPParametricAtlas::parentPoint( const topology::Vertex& v ) const
    {
        // We do some finagling to get the appropriate points from the underlying atlases.

        const auto [source_dart, line_dart, tp_pos] = mMap->unflatten( v.dart() );

        const ParentPoint line_ppt = [&]() {
            if( std::to_underlying( tp_pos ) >= mMap->dim() )
                return ParentPoint( simplexDomain( 1 ), Eigen::Matrix<double, 1, 1>( 1.0 ), { true, false } );
            else
                return ParentPoint( simplexDomain( 1 ), Eigen::Matrix<double, 1, 1>( 0.0 ), { false, true } );
        }();

        const ParentPoint source_ppt = [&]() {
            if( mMap->dim() == 3 )
            {
                switch( tp_pos )
                {
                    case topology::TPCombinatorialMap::TPDartPos::DartPos0:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos2:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos3:
                        return mSourceParam->parentPoint( source_dart );
                    default: return mSourceParam->parentPoint( phi( mSourceParam->cmap(), 1, source_dart ).value() );
                }
            }
            else
            {
                switch( tp_pos )
                {
                    case topology::TPCombinatorialMap::TPDartPos::DartPos0:
                    case topology::TPCombinatorialMap::TPDartPos::DartPos3:
                        return ParentPoint( simplexDomain( 1 ), Eigen::Matrix<double, 1, 1>( 0.0 ), { false, true } );
                    default:
                        return ParentPoint( simplexDomain( 1 ), Eigen::Matrix<double, 1, 1>( 1.0 ), { true, false } );
                }
            }
        }();

        return tensorProduct( source_ppt, line_ppt );
    }

    Vector6dMax TPParametricAtlas::parametricLengths( const topology::Cell& c ) const
    {
        const auto [source_dart, line_dart, tp_pos] = mMap->unflatten( c.dart() );
        const Vector6dMax source_lengths = mSourceParam->parametricLengths( topology::Cell( source_dart, mSourceParam->cmap().dim() ) );
        const Vector6dMax line_lengths = mLineParam->parametricLengths( topology::Cell( line_dart, mLineParam->cmap().dim() ) );
        return ( Vector6dMax( source_lengths.size() + line_lengths.size() ) << source_lengths, line_lengths ).finished();
    }

    SmallVector<topology::Cell, 12> cornerCells( const TPParametricAtlas& atlas, const uint cell_dim )
    {
        using namespace topology;
        const TPCombinatorialMap& cmap = atlas.cmap();
        const auto comps = tensorProductComponentCMaps( cmap );
        if( comps.size() != cmap.dim() )
            throw std::runtime_error( "cornerCells only works on cube-like TP regions" );

        if( cell_dim >= cmap.dim() )
            throw std::runtime_error( "cornerCells requires a cell dim less than the cmap dim" );

        const util::IndexVec lengths = [&comps]() {
            util::IndexVec out;
            for( const auto& comp : comps ) out.push_back( cellCount( *comp, comp->dim() ) );
            return out;
        }();

        const auto iterate_combinations =
            []( const size_t dim, const size_t num, const std::function<void( const util::IndexVec& )>& callback ) {
                SmallVector<bool, 3> a( num, true );
                a.resize( dim, false );
                do
                {
                    util::IndexVec out;
                    for( size_t i = 0; i < dim; i++ )
                        if( a.at( i ) ) out.push_back( i );
                    callback( out );
                } while( std::prev_permutation( a.begin(), a.end() ) );
            };

        SmallVector<Cell, 12> out;
        for( size_t i = 0, n = pow( 2, lengths.size() - cell_dim ); i < n; i++ )
        {
            const util::IndexVec partial_bdry = util::unflatten( i, util::IndexVec( lengths.size() - cell_dim, 2 ) );

            iterate_combinations( cmap.dim(), cell_dim, [&]( const util::IndexVec& zeros_to_add ) {
                util::IndexVec which_bdry = partial_bdry;
                for( const size_t coord : zeros_to_add )
                    which_bdry.insert( std::next( which_bdry.begin(), coord ), 0 );

                SmallVector<Dart, 3> unflat_darts;
                std::transform( which_bdry.begin(),
                                which_bdry.end(),
                                lengths.begin(),
                                std::back_inserter( unflat_darts ),
                                []( const size_t which_bdry, const size_t len ) { return Dart( which_bdry * ( len - 1 ) ); } );

                const Dart flat_dart = flattenFull( cmap, unflat_darts );

                const BaryCoordIsZeroVec bdry_zero_vec = [&]() {
                    BaryCoordIsZeroVec out;
                    for( size_t i = 0, partial_i = 0; i < cmap.dim(); i++ )
                    {
                        if( std::ranges::find( zeros_to_add, i ) == zeros_to_add.end() )
                        {
                            const bool is_zero = partial_bdry.at( partial_i ) == 0;
                            out.push_back( not is_zero );
                            out.push_back( is_zero );
                            partial_i++;
                        }
                        else
                        {
                            out.push_back( false );
                            out.push_back( false );
                        }
                    }
                    return out;
                }();

                const bool found_one = not iterateDartsOfCell( cmap, Cell( flat_dart, cmap.dim() ), [&]( const Dart& d ) {
                    const Cell c( d, cell_dim );
                    if( parentDomainBoundary( atlas, c ) == bdry_zero_vec )
                    {
                        out.push_back( c );
                        return false;
                    }
                    return true;
                } );

                if( not found_one ) std::runtime_error( "Issue here!" );
            } );
        }
        return out;
    }
}