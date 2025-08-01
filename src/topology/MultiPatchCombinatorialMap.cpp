#include <MultiPatchCombinatorialMap.hpp>
#include <IndexOperations.hpp>
#include <TPCombinatorialMap.hpp>
#include <utility>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>
#include <algorithm>
#include <numeric>

#include <Logging.hpp>

namespace topology
{
    using TPPermutation = MultiPatchCombinatorialMap::TPPermutation;
    using TPDartPos = TPCombinatorialMap::TPDartPos;

    static constexpr std::array<std::array<Dart::IndexType, 4>, 6> local_darts_from_bdry_ids_3d{
        std::array<Dart::IndexType, 4>{ 7, 8, 9, 10 },
        std::array<Dart::IndexType, 4>{ 19, 20, 21, 22 },
        std::array<Dart::IndexType, 4>{ 13, 14, 15, 16 },
        std::array<Dart::IndexType, 4>{ 1, 2, 3, 4 },
        std::array<Dart::IndexType, 4>{ 5, 23, 17, 11 },
        std::array<Dart::IndexType, 4>{ 0, 6, 12, 18 } };

    // Boundaries are specified by the entry in BoundaryBoolVec that would be true.
    // This results in the ordering S1, S0, T1, T0, ...
    static constexpr std::array<size_t, 4> local_dart_to_bdry_2d{ 3, 0, 2, 1 };
    static constexpr std::array<Dart::IndexType, 4> local_dart_from_bdry_2d{ 1, 3, 2, 0 };

    Dart flipDart( const util::IndexVec& tp_lengths, const SmallVector<Dart, 2> from_darts, const size_t id )
    {
        return Dart( tp_lengths.at( id ) - from_darts.at( id ).id() - 1 );
    }

    SmallVector<Dart, 2> flipNecessaryDarts( const util::IndexVec& bdry_lengths,
                                             const SmallVector<Dart, 2>& bdry_unflat,
                                             const size_t which_boundary )
    {
        if( bdry_unflat.size() == 1 )
        {
            if( which_boundary == 1 or which_boundary == 2 )
                return { flipDart( bdry_lengths, bdry_unflat, 0 ) };
            else
                return bdry_unflat;
        }
        else
        {
            switch( which_boundary )
            {
                case 0:
                case 3:
                case 4: return { flipDart( bdry_lengths, bdry_unflat, 0 ), bdry_unflat.at( 1 ) };
                case 1:
                case 2:
                case 5: return bdry_unflat;
                default: throw std::runtime_error( "Bad boundary" );
            }
        }
    };

    std::pair<SmallVector<Dart, 2>, TPDartPos> permuteAcrossBoundary( const SmallVector<Dart, 2>& from_darts,
                                                                      const TPDartPos& from_pos,
                                                                      const TPPermutation& permute,
                                                                      const util::IndexVec& tp_lengths )
    {
        if( tp_lengths.size() != from_darts.size() )
            throw std::runtime_error( "tp_lengths size must match from_darts size" );

        if( from_darts.size() != ( permute == TPPermutation::Flip1d ? 1 : 2 ) )
            throw std::runtime_error( "Bad inputs to permuteAcrossBoundary." );

        const TPDartPos permuted_pos =
            TPDartPos( ( std::to_underlying( permute ) + ( 4 - std::to_underlying( from_pos ) ) ) % 4 );

        switch( permute )
        {
            case TPPermutation::Flip1d:
                // 1d boundary, just flip the dart id
                return { { flipDart( tp_lengths, from_darts, 0 ) }, from_pos };
            case TPPermutation::ZeroToZero:
                return { { flipDart( tp_lengths, from_darts, 0 ), from_darts.at( 1 ) }, permuted_pos };
            case TPPermutation::ZeroToOne:
                return { { flipDart( tp_lengths, from_darts, 1 ), flipDart( tp_lengths, from_darts, 0 ) },
                         permuted_pos };
            case TPPermutation::ZeroToTwo:
                return { { from_darts.at( 0 ), flipDart( tp_lengths, from_darts, 1 ) }, permuted_pos };
            case TPPermutation::ZeroToThree: return { { from_darts.at( 1 ), from_darts.at( 0 ) }, permuted_pos };
            default: return { {}, TPDartPos::DartPos0 };// Unreachable; to silence warning on ubuntu
        }
    }

    size_t whichBoundary( const TPCombinatorialMap& cmap, const Dart& d )
    {
        if( cmap.dim() == 3 )
        {
            const Dart::IndexType hex_local_d = d.id() % 24;
            for( size_t i = 0; i < local_darts_from_bdry_ids_3d.size(); i++ )
            {
                if( std::ranges::find( local_darts_from_bdry_ids_3d.at( i ), hex_local_d ) !=
                    local_darts_from_bdry_ids_3d.at( i ).end() )
                {
                    return i;
                }
            }
            throw std::runtime_error( "Dart on no boundary!" );
        }
        else
            return local_dart_to_bdry_2d.at( d.id() % 4 );
    }

    std::pair<SmallVector<Dart, 2>, TPDartPos> dropToBoundary( const TPCombinatorialMap& cmap,
                                                               const Dart& d,
                                                               const size_t which_boundary,
                                                               const util::IndexVec& bdry_lengths )
    {
        const TPDartPos dart_pos = [&]() {
            if( cmap.dim() == 3 )
            {
                const auto& local_darts_from_bdry = local_darts_from_bdry_ids_3d.at( which_boundary );
                const auto it = std::find( local_darts_from_bdry.begin(), local_darts_from_bdry.end(), d.id() % 24 );
                if( it == local_darts_from_bdry.end() ) throw std::runtime_error( "Dart not in provided boundary!" );
                return TPDartPos( std::distance( local_darts_from_bdry.begin(), it ) );
            }
            else
                return TPDartPos( local_dart_to_bdry_2d.at( d.id() % 4 ) );
        }();

        const SmallVector<Dart, 3> full_unflat = unflattenFull( cmap, d ).unflat_darts;

        // Drop the correct dart to get to the boundary
        SmallVector<Dart, 2> bdry_unflat;
        for( size_t i = 0; i < full_unflat.size(); i++ )
        {
            if( i != which_boundary / 2 ) bdry_unflat.push_back( full_unflat.at( i ) );
        }

        return { flipNecessaryDarts( bdry_lengths, bdry_unflat, which_boundary ), dart_pos };
    }

    util::IndexVec boundaryLengths( const TPCombinatorialMap& cmap, const size_t which_boundary )
    {
        const auto comps = tensorProductComponentCMaps( cmap );

        util::IndexVec out;
        for( size_t i = 0; i < comps.size(); i++ )
        {
            if( i == which_boundary / 2 ) continue;
            const auto& comp = comps.at( i );
            out.push_back( cellCount( *comp, comp->dim() ) );
        }

        return out;
    }

    size_t perpLength( const TPCombinatorialMap& cmap, const size_t which_boundary )
    {
        const auto comps = tensorProductComponentCMaps( cmap );

        for( size_t i = 0; i < comps.size(); i++ )
        {
            const auto& comp = comps.at( i );
            if( i == which_boundary / 2 ) return cellCount( *comp, comp->dim() );
        }

        throw std::runtime_error( "Bad boundary specified" );
    }

    Dart raiseFromBoundary( const TPCombinatorialMap& cmap,
                            const SmallVector<Dart, 2>& from_darts,
                            const TPDartPos& from_pos,
                            const size_t which_boundary,
                            const util::IndexVec& bdry_lengths )
    {
        const auto bdry_darts = flipNecessaryDarts( bdry_lengths, from_darts, which_boundary );
        SmallVector<Dart, 3> unflat_darts;
        for( size_t i = 0; i <= bdry_darts.size(); i++ )
        {
            if( i == which_boundary / 2 )
            {
                unflat_darts.push_back( [&]() {
                    if( which_boundary % 2 == 1 )
                        return Dart( 0 );
                    else
                        return Dart( perpLength( cmap, which_boundary ) - 1 );
                }() );
            }
            if( i < bdry_darts.size() ) unflat_darts.push_back( bdry_darts.at( i ) );
        }

        const Dart::IndexType local_dart_id =
            cmap.dim() == 2 ? local_dart_from_bdry_2d.at( which_boundary )
                            : local_darts_from_bdry_ids_3d.at( which_boundary ).at( std::to_underlying( from_pos ) );

        return Dart( flattenFull( cmap, unflat_darts ).id() + local_dart_id );
    }

    using ConstituentSide = MultiPatchCombinatorialMap::ConstituentSide;
    std::map<ConstituentSide, std::pair<TPPermutation, ConstituentSide>> initializeInterMapConnections(
        const std::vector<std::shared_ptr<const TPCombinatorialMap>>& constituents,
        const std::map<std::pair<size_t, Dart>, std::pair<size_t, Dart>>& connections )
    {
        if( constituents.size() == 0 ) throw std::runtime_error( "Cannot create multipatch of zero patches" );

        std::map<ConstituentSide, std::pair<TPPermutation, ConstituentSide>> out;
        for( const auto& connection : connections )
        {
            const TPCombinatorialMap& left_cmap = *constituents.at( connection.first.first );
            const TPCombinatorialMap& right_cmap = *constituents.at( connection.second.first );
            const size_t left_bdry = whichBoundary( left_cmap, connection.first.second );
            const size_t right_bdry = whichBoundary( right_cmap, connection.second.second );

            const ConstituentSide left_side{ connection.first.first, left_bdry };
            const ConstituentSide right_side{ connection.second.first, right_bdry };

            const TPPermutation perm = [&]() -> TPPermutation {
                if( constituents.front()->dim() == 2 ) return TPPermutation::Flip1d;

                const TPDartPos left_pos =
                    dropToBoundary(
                        left_cmap, connection.first.second, left_bdry, boundaryLengths( left_cmap, left_bdry ) )
                        .second;
                const TPDartPos right_pos =
                    dropToBoundary(
                        right_cmap, connection.second.second, right_bdry, boundaryLengths( right_cmap, right_bdry ) )
                        .second;

                return TPPermutation( ( std::to_underlying( left_pos ) + std::to_underlying( right_pos ) ) % 4 );
            }();

            out.emplace( left_side, std::pair<TPPermutation, ConstituentSide>{ perm, right_side } );
            out.emplace( right_side, std::pair<TPPermutation, ConstituentSide>{ perm, left_side } );
        }

        return out;
    }

    MultiPatchCombinatorialMap::MultiPatchCombinatorialMap(
        const std::vector<std::shared_ptr<const TPCombinatorialMap>>& constituents,
        const std::map<std::pair<size_t, Dart>, std::pair<size_t, Dart>>& connections )
        : MultiPatchCombinatorialMap( constituents, initializeInterMapConnections( constituents, connections ) )
    {}

    MultiPatchCombinatorialMap::MultiPatchCombinatorialMap(
        const std::vector<std::shared_ptr<const TPCombinatorialMap>>& constituents,
        const std::map<ConstituentSide, std::pair<TPPermutation, ConstituentSide>>& connections )
        : mSubMaps( constituents ), mInterMapConnections( connections ), mRanges( initializeRanges( constituents ) )
    {}

    std::pair<size_t, Dart> MultiPatchCombinatorialMap::toLocalDart( const Dart& global_dart ) const
    {
        return mRanges.toLocalDart( global_dart );
    }

    Dart MultiPatchCombinatorialMap::toGlobalDart( const size_t patch_id, const Dart& local_dart ) const
    {
        return mRanges.toGlobalDart( patch_id, local_dart );
    }

    std::optional<Dart> MultiPatchCombinatorialMap::phi( const int i, const Dart& d ) const
    {
        const auto [patch_id, local_d] = toLocalDart( d );

        const TPCombinatorialMap& patch_cmap = *mSubMaps.at( patch_id );

        return topology::phi( patch_cmap, i, local_d )
            .transform( [&]( const Dart& local_phi ) { return toGlobalDart( patch_id, local_phi ); } )
            .or_else( [&]() -> std::optional<Dart> {
                const size_t bdry_id = whichBoundary( patch_cmap, local_d );
                const auto it = mInterMapConnections.find( { patch_id, bdry_id } );
                if( it == mInterMapConnections.end() ) return {};

                const std::pair<TPPermutation, ConstituentSide>& connection = it->second;
                const util::IndexVec bdry_lengths = boundaryLengths( patch_cmap, bdry_id );

                const std::pair<SmallVector<Dart, 2>, TPDartPos> bdry_darts =
                    dropToBoundary( patch_cmap, local_d, bdry_id, bdry_lengths );
                const std::pair<SmallVector<Dart, 2>, TPDartPos> next_bdry_darts =
                    permuteAcrossBoundary( bdry_darts.first, bdry_darts.second, connection.first, bdry_lengths );

                const Dart next_local_d = raiseFromBoundary(
                    *mSubMaps.at( connection.second.constituent_id ),
                    next_bdry_darts.first,
                    next_bdry_darts.second,
                    connection.second.side_id,
                    boundaryLengths( *mSubMaps.at( connection.second.constituent_id ), connection.second.side_id ) );

                return toGlobalDart( connection.second.constituent_id, next_local_d );
            } );
    }

    Dart::IndexType MultiPatchCombinatorialMap::maxDartId() const
    {
        return mRanges.maxDartId();
    }

    uint MultiPatchCombinatorialMap::dim() const
    {
        return mSubMaps.front()->dim();
    }

    bool MultiPatchCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
    {
        for( size_t patch_id = 0; patch_id < mSubMaps.size(); patch_id++ )
        {
            const bool continue_iter =
                topology::iterateDartsWhile( *mSubMaps.at( patch_id ), [&]( const Dart& local_dart ) {
                    return callback( toGlobalDart( patch_id, local_dart ) );
                } );
            if( not continue_iter ) return false;
        }
        return true;
    }

    bool MultiPatchCombinatorialMap::iterateCellsWhile( const uint cell_dim,
                                                        const std::function<bool( const Cell& )>& callback ) const
    {
        if( cell_dim == dim() )
        {
            for( size_t patch_id = 0; patch_id < mSubMaps.size(); patch_id++ )
            {
                const bool continue_iter =
                    topology::iterateCellsWhile( *mSubMaps.at( patch_id ), cell_dim, [&]( const Cell& local_cell ) {
                        return callback( Cell( toGlobalDart( patch_id, local_cell.dart() ), local_cell.dim() ) );
                    } );

                if( not continue_iter ) return false;
            }
            return true;
        }
        else
        {
            GlobalCellMarker m( *this, cell_dim );
            for( size_t patch_id = 0; patch_id < mSubMaps.size(); patch_id++ )
            {
                const bool continue_iter =
                    topology::iterateCellsWhile( *mSubMaps.at( patch_id ), cell_dim, [&]( const Cell& local_cell ) {
                        const Cell global_cell( toGlobalDart( patch_id, local_cell.dart() ), local_cell.dim() );
                        if( not m.isMarked( global_cell ) )
                        {
                            m.mark( *this, global_cell );
                            return callback( global_cell );
                        }
                        return true;
                    } );

                if( not continue_iter ) return false;
            }
            return true;
        }
    }

    std::optional<IndexingFunc> MultiPatchCombinatorialMap::indexing( const uint ) const
    {
        return std::nullopt;
        // std::vector<IndexingFunc> indexing_funcs;
        // indexing_funcs.reserve( mSubMaps.size() );
        // for( const auto& sub : mSubMaps )
        // {
        //     const std::optional<IndexingFunc> sub_idx = sub->indexing( cell_dim );
        //     if( not sub_idx.has_value() ) return std::nullopt;
        //     indexing_funcs.push_back( sub_idx.value() );
        // }
        // return [&this, indexing_funcs]( const Cell& c ) {
        //     const Cell local_cell( toLocalDart)
        // }
    }

    std::optional<size_t> MultiPatchCombinatorialMap::cellCount( const uint cell_dim ) const
    {
        if( cell_dim != dim() ) return std::nullopt;

        return std::transform_reduce( mSubMaps.begin(),
                                      mSubMaps.end(),
                                      0,
                                      std::plus<>(),
                                      [&cell_dim]( const std::shared_ptr<const TPCombinatorialMap>& m ) {
                                          return topology::cellCount( *m, cell_dim );
                                      } );
    }

    bool MultiPatchCombinatorialMap::ConstituentSide::operator<(
        const MultiPatchCombinatorialMap::ConstituentSide& o ) const
    {
        if( constituent_id == o.constituent_id )
            return side_id < o.side_id;
        else
            return constituent_id < o.constituent_id;
    }

    std::ostream& operator<<( std::ostream& o, const MultiPatchCombinatorialMap::ConstituentSide& cs )
    {
        o << "{" << cs.constituent_id << ", " << cs.side_id << "}";
        return o;
    }
    std::ostream& operator<<( std::ostream& o, const MultiPatchCombinatorialMap::TPPermutation& tpp )
    {
        switch( tpp )
        {
            case MultiPatchCombinatorialMap::TPPermutation::Flip1d:
                o << "Flip1d";
                break;
            case MultiPatchCombinatorialMap::TPPermutation::ZeroToOne:
                o << "ZeroToOne";
                break;
            case MultiPatchCombinatorialMap::TPPermutation::ZeroToTwo:
                o << "ZeroToTwo";
                break;
            case MultiPatchCombinatorialMap::TPPermutation::ZeroToThree:
                o << "ZeroToThree";
                break;
            case MultiPatchCombinatorialMap::TPPermutation::ZeroToZero:
                o << "ZeroToZero";
                break;
        }
        return o;
    }

    MultiPatchCombinatorialMap::InternalConnectionsMap
        connectionsOfSweptMultipatch( const MultiPatchCombinatorialMap::InternalConnectionsMap& connections_2d )
    {
        using InternalConnectionsMap = MultiPatchCombinatorialMap::InternalConnectionsMap;
        using ConstituentSide = MultiPatchCombinatorialMap::ConstituentSide;
        using TPPermutation = MultiPatchCombinatorialMap::TPPermutation;

        InternalConnectionsMap out;
        for( const auto& [left_side, perm_and_right_side] : connections_2d )
        {
            if( perm_and_right_side.first != TPPermutation::Flip1d )
                throw std::runtime_error( "Can only sweep multipatch connections of 2d maps." );

            // All constituents have the same orientation in this situation, so use ZeroToZero.
            out.emplace( left_side, std::pair<TPPermutation, ConstituentSide>{ TPPermutation::ZeroToZero, perm_and_right_side.second } );
        }

        return out;
    }

    MultiPatchCombinatorialMap blockLayout( const MultiPatchCombinatorialMap& cmap )
    {
        const std::shared_ptr<const TPCombinatorialMap> single_elem_patch = [&]() {
            const auto cmap1d = std::make_shared<const CombinatorialMap1d>( 1 );
            const auto cmap2d = std::make_shared<const TPCombinatorialMap>( cmap1d, cmap1d );
            if( cmap.dim() == 3 ) {
                return std::make_shared<const TPCombinatorialMap>( cmap2d, cmap1d );
            } else {
                return cmap2d;
            }
        }();

        const MultiPatchCombinatorialMap block_layout(
            std::vector<std::shared_ptr<const TPCombinatorialMap>>( cmap.constituents().size(), single_elem_patch ),
            cmap.connections() );
        return block_layout;
    }
} // namespace topology