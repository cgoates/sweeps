#include <HierarchicalTPCombinatorialMap.hpp>
#include <GlobalCellMarker.hpp>
#include <IndexOperations.hpp>
#include <ranges>
#include <numeric>

using namespace topology;

namespace topology
{
bool checkForNoAncestor( const FullyUnflattenedDart& unflat, const size_t n_darts_per_ancestor )
{
    using TPDartPos = TPCombinatorialMap::TPDartPos;

    const auto check_dart = [&]( const size_t idx, const bool add_one = false ) {
        return ( unflat.unflat_darts.at( idx ).id() + ( add_one ? 1 : 0 ) ) % n_darts_per_ancestor != 0;
    };
    constexpr bool add_one = true;
    if( unflat.unflat_darts.size() == 2 ) // 2D case
    {
        switch( unflat.dart_pos.at( 0 ) )
        {
            case TPCombinatorialMap::TPDartPos::DartPos0:
                return check_dart( 1 );
            case TPCombinatorialMap::TPDartPos::DartPos1:
                return check_dart( 0, add_one );
            case TPCombinatorialMap::TPDartPos::DartPos2:
                return check_dart( 1, add_one );
            case TPCombinatorialMap::TPDartPos::DartPos3:
                return check_dart( 0 );
            default:
                throw std::runtime_error( "Invalid TPDartPos for 2d" );
        }
    }
    else if( unflat.unflat_darts.size() == 3 ) // 3D case
    {
        switch( unflat.dart_pos.at( 0 ) )
        {
            case TPDartPos::DartPos0:
                switch( unflat.dart_pos.at( 1 ) )
                {
                    case TPDartPos::DartPos0:
                    case TPDartPos::DartPos1:
                        return check_dart( 1 ) or check_dart( 2 );
                    case TPDartPos::DartPos2:
                        return check_dart( 0 ) or check_dart( 1 );
                    case TPDartPos::DartPos4:
                        return check_dart( 0, add_one ) or check_dart( 1 );
                    case TPDartPos::DartPos3:
                    case TPDartPos::DartPos5:
                        return check_dart( 1 ) or check_dart( 2, add_one );
                    default:
                        throw std::runtime_error( "Invalid TPDartPos for 3d" );
                }
            case TPDartPos::DartPos1:
                switch( unflat.dart_pos.at( 1 ) )
                {
                    case TPDartPos::DartPos0:
                    case TPDartPos::DartPos1:
                        return check_dart( 0, add_one ) or check_dart( 2 );
                    case TPDartPos::DartPos2:
                        return check_dart( 0, add_one ) or check_dart( 1 );
                    case TPDartPos::DartPos4:
                        return check_dart( 0, add_one ) or check_dart( 1, add_one );
                    case TPDartPos::DartPos3:
                    case TPDartPos::DartPos5:
                        return check_dart( 0, add_one ) or check_dart( 2, add_one );
                    default:
                        throw std::runtime_error( "Invalid TPDartPos for 3d" );
                }
            case TPDartPos::DartPos2:
                switch( unflat.dart_pos.at( 1 ) )
                {
                    case TPDartPos::DartPos0:
                    case TPDartPos::DartPos1:
                        return check_dart( 1, add_one ) or check_dart( 2 );
                    case TPDartPos::DartPos2:
                        return check_dart( 0, add_one ) or check_dart( 1, add_one );
                    case TPDartPos::DartPos4:
                        return check_dart( 0 ) or check_dart( 1, add_one );
                    case TPDartPos::DartPos3:
                    case TPDartPos::DartPos5:
                        return check_dart( 1, add_one ) or check_dart( 2, add_one );
                    default:
                        throw std::runtime_error( "Invalid TPDartPos for 3d" );
                }
            case TPDartPos::DartPos3:
                switch( unflat.dart_pos.at( 1 ) )
                {
                    case TPDartPos::DartPos0:
                    case TPDartPos::DartPos1:
                        return check_dart( 0 ) or check_dart( 2 );
                    case TPDartPos::DartPos2:
                        return check_dart( 0 ) or check_dart( 1, add_one );
                    case TPDartPos::DartPos4:
                        return check_dart( 0 ) or check_dart( 1 );
                    case TPDartPos::DartPos3:
                    case TPDartPos::DartPos5:
                        return check_dart( 0 ) or check_dart( 2, add_one );
                    default:
                        throw std::runtime_error( "Invalid TPDartPos for 3d" );
                }

            default:
                throw std::runtime_error( "Invalid TPDartPos for 2d" );
        }
    }
    throw std::invalid_argument( "Bad unflattened dart passed to checkForNoAncestor" );
}

bool checkForNoAncestor( const TPCombinatorialMap& tp_map, const Dart& d, const size_t n_darts_per_ancestor )
{
    if( n_darts_per_ancestor == 0 )
        return true;
    return checkForNoAncestor( unflattenFull( tp_map, d ), n_darts_per_ancestor );
}
}

HierarchicalTPCombinatorialMap::HierarchicalTPCombinatorialMap(
    const std::vector<std::shared_ptr<const TPCombinatorialMap>>& refinement_levels,
    const std::vector<std::vector<Cell>>& leaf_elements )
    : mRefinementLevels( refinement_levels ),
      mRanges( initializeRanges( refinement_levels ) ),
      mLeafDarts( mRanges.maxDartId() + 1, false ),
      mUnrefinedDarts( mRanges.maxDartId() + 1, false )
{
    const auto get_tp_lengths = []( const TPCombinatorialMap& cmap ) -> util::IndexVec {
        const auto cmaps = tensorProductComponentCMaps( cmap );
        util::IndexVec out;
        std::transform( cmaps.begin(), cmaps.end(), std::back_inserter( out ), []( const auto& c ) { return topology::cellCount( *c, 1 ); } );
        return out;
    };

    // For now I am assuming uniform refinement only.
    util::IndexVec lower_sizes = get_tp_lengths( *refinement_levels.at( 0 ) );
    for( size_t level = 1; level < refinement_levels.size(); level++ )
    {
        const util::IndexVec higher_sizes = get_tp_lengths( *refinement_levels.at( level ) );
        mRefinementRatios.push_back( higher_sizes.front() / lower_sizes.front() );

        if( higher_sizes.front() % lower_sizes.front() != 0 )
            throw std::invalid_argument( "Refinement levels must be uniform refinements of each other." );

        for( size_t i = 0; i < lower_sizes.size(); i++ )
        {
            if( higher_sizes.at( i ) / lower_sizes.at( i ) != mRefinementRatios.back() or higher_sizes.at( i ) % lower_sizes.at( i ) != 0 )
                throw std::invalid_argument( "Refinement levels must be uniform refinements of each other." );
        }
        lower_sizes = higher_sizes;
    }

    std::set<Dart> not_leaf;

    const auto mark_as_leaf = [&]( const size_t level, const Dart& level_d ) {
        const Dart global_d = mRanges.toGlobalDart( level, level_d );
        mLeafDarts.at( global_d.id() ) = true;

        // Unmark ancestors that are marked as leaf.
        iterateAncestors( global_d, [&]( const Dart& ancestor_global_d ) {
            not_leaf.emplace( ancestor_global_d );
            return true;
        } );
    };

    const auto has_leaf_ancestor = [&]( const size_t level_ii, const Dart& level_d ) {
        const Dart global_d = mRanges.toGlobalDart( level_ii, level_d );
        return not iterateAncestors( global_d, [&]( const Dart& ancestor_global_d ) {
            return not mLeafDarts.at( ancestor_global_d.id() );
        } );
    };

    // Mark unrefined and leaf darts
    for( size_t level = 0; level < refinement_levels.size(); level++ )
    {
        for( const auto& leaf_elem : leaf_elements.at( level ) )
        {
            iterateDartsOfCell( *refinement_levels.at( level ), leaf_elem, [&]( const Dart& level_d ) {
                mUnrefinedDarts.at( mRanges.toGlobalDart( level, level_d ).id() ) = true;
                const auto maybe_phi = topology::phi( *refinement_levels.at( level ), dim(), level_d );

                mark_as_leaf( level, level_d );
                if( maybe_phi.has_value() )
                    mark_as_leaf( level, maybe_phi.value() );
                if( dim() == 3 )
                {
                    iterateDartsOfCell( *refinement_levels.at( level ), Edge( level_d ), [&]( const Dart& other_level_d ) {
                        if( has_leaf_ancestor( level, other_level_d ) )
                        {
                            mark_as_leaf( level, other_level_d );
                            const auto maybe_other_phi = topology::phi( *refinement_levels.at( level ), dim(), other_level_d );
                            if( maybe_other_phi.has_value() )
                                mark_as_leaf( level, other_level_d );
                        }
                        return true;
                    } );
                }

                return true;
            } );
        }
    }

    // Unmark ancestors that are marked as leaf.
    for( const auto& d : not_leaf )
    {
        mLeafDarts.at( d.id() ) = false;
    }

    // Build phi1 and phi-1 ops
    for( size_t level = 0; level < refinement_levels.size(); level++ )
    {
        const TPCombinatorialMap& level_cmap = *refinement_levels.at( level );
        for( const auto& leaf_elem : leaf_elements.at( level ) )
        {
            const auto phi_chain = [&level_cmap, &level, this]( const Face& leaf_face ) {
                const bool all_darts_are_leaves = iterateDartsOfRestrictedCell( level_cmap, leaf_face, 3, [&]( const Dart& d ) {
                    const Dart global_d = mRanges.toGlobalDart( level, d );
                    return mLeafDarts.at( global_d.id() );
                } );

                if( not all_darts_are_leaves )
                {
                    std::optional<Dart> first_leaf;
                    std::optional<Dart> previous_leaf;
                    bool perform_phi1_chain = true;
                    const Dart global_leaf_face_dart = mRanges.toGlobalDart( level, leaf_face.dart() );
                    Dart d = global_leaf_face_dart;
                    std::map<Dart, Dart> face_local_phi1s;
                    do
                    {
                        std::optional<Dart> has_leaf_phi1;
                        iterateLeafDescendants( d, [&]( const Dart& descendant ) {
                            if( previous_leaf )
                            {
                                face_local_phi1s.insert( { previous_leaf.value(), descendant } );
                            }
                            previous_leaf.emplace( descendant );
                            if( not first_leaf ) first_leaf.emplace( descendant );

                            if( not has_leaf_phi1 )
                            {
                                const auto [desc_level, desc_level_d] = mRanges.toLocalDart( descendant );
                                const Dart phi1 = topology::phi( *mRefinementLevels.at( desc_level ), 1, desc_level_d ).value();
                                if( mLeafDarts.at( mRanges.toGlobalDart( desc_level, phi1 ).id() ) )
                                {
                                    has_leaf_phi1.emplace( descendant );
                                }
                            }

                            return true;
                        } );

                        if( has_leaf_phi1.has_value() and previous_leaf.has_value() and
                            previous_leaf.value() != has_leaf_phi1.value() )
                        {
                            // We shouldn't do the phi chain in this case.
                            // This means that the face is actually multiple faces.
                            // *------------------------*
                            // |                        |
                            // |                        |
                            // |                        |
                            // |                        |
                            // |          |             |
                            // |        b |             |
                            // |      a   |       c     |
                            // | *------- *  *--------  |
                            // *------------------------*
                            // a has a leaf phi1 in its level, means that b is a leaf dart, and therefore
                            // the face is actually refined.  If it were c that had a leaf phi1, then
                            // that doesn't actually tell us that the face is refined.
                            // If the face is refined, Phi 1 and -1 operations will be handled by a higher
                            // refinement level, as there are higher level cells on the phi3.
                            perform_phi1_chain = false;
                            break;
                        }
                        d = mRanges.toGlobalDart( level, topology::phi( level_cmap, 1, mRanges.toLocalDart( d ).second ).value() );
                    } while ( d != global_leaf_face_dart );

                    if( perform_phi1_chain )
                    {
                        if( first_leaf )
                        {
                            face_local_phi1s.insert( { previous_leaf.value(), first_leaf.value() } );
                        }

                        for( const auto& [key, value] : face_local_phi1s )
                        {
                            mPhiOnes.insert( { key, value } );
                            mPhiMinusOnes.insert( { value, key } );
                        }
                    }
                }
            };

            // In 3d we need to do this for the faces of the leaf elements, in 2d, the elements themselves.
            iterateAdjacentCells( level_cmap, leaf_elem, 2, [&]( const topology::Face& leaf_face ) {
                phi_chain( leaf_face );
                if( dim() == 3 )
                {
                    const auto maybe_phi3 = topology::phi( level_cmap, 3, leaf_face.dart() );
                    if( maybe_phi3.has_value() )
                    {
                        phi_chain( maybe_phi3.value() );
                    }
                }
                return true;
            } );
        }
    }
}

std::optional<Dart> HierarchicalTPCombinatorialMap::phi( const int i, const Dart& d ) const
{
    const auto [level, level_d ] = mRanges.toLocalDart( d );
    return topology::phi( *mRefinementLevels.at( level ), i, level_d ).and_then( [&]( const Dart& level_phi ) -> std::optional<Dart> {
        const Dart level_global_phi = mRanges.toGlobalDart( level, level_phi );

        if( mLeafDarts.at( level_global_phi.id() ) ) return level_global_phi;

        if( i == 2 and dim() == 3 )
        {
            // Sometimes a phi2 in 3d needs to be a phi(2,3,2)
            return topology::phi( *mRefinementLevels.at( level ), { 2, 3, 2 }, level_d ).and_then( [&]( const Dart& phi_323 ) -> std::optional<Dart> {
                return mRanges.toGlobalDart( level, phi_323 );
            } );
        }
        else
        {
            if( i > 1 ) throw std::runtime_error( "Missing phi(dim) operation!" );
            const std::map<Dart, Dart>& phi_map = i == 1 ? mPhiOnes : mPhiMinusOnes;
            const auto it = phi_map.find( d );
            if( it != phi_map.end() ) return it->second;
            else throw std::runtime_error( "Phi " + std::to_string( i ) + " operation unfound!" );
        }
    } );
}

Dart::IndexType HierarchicalTPCombinatorialMap::maxDartId() const
{
    return mRanges.maxDartId();
}

uint HierarchicalTPCombinatorialMap::dim() const
{
    return mRefinementLevels.front()->dim();
}

bool HierarchicalTPCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    for( size_t i = 0; i < mLeafDarts.size(); i++ )
    {
        if( mLeafDarts.at( i ) )
        {
            if( not callback( Dart( i ) ) ) return false;
        }
    }
    return true;
}

bool HierarchicalTPCombinatorialMap::iterateCellsWhile( const uint cell_dim,
                                                        const std::function<bool( const Cell& )>& callback ) const
{
    // TODO: Do more here for better performance?
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

std::optional<IndexingFunc> HierarchicalTPCombinatorialMap::indexing( const uint ) const
{
    return std::nullopt;
}

std::optional<size_t> HierarchicalTPCombinatorialMap::cellCount( const uint ) const
{
    return std::nullopt; // TODO: Do more here for better performance?
}

std::pair<size_t, Dart> HierarchicalTPCombinatorialMap::unrefinedAncestorDartOfCell( const Cell& leaf_c ) const
{
    const Dart& leaf_d = leaf_c.dart();
    if( not mLeafDarts.at( leaf_d.id() ) )
        throw std::runtime_error( "unrefinedAncestorDartOfCell requires a leaf dart as input" );

    if( mUnrefinedDarts.at( leaf_d.id() ) ) return mRanges.toLocalDart( leaf_d );

    std::optional<Dart> global_ancestor_d;

    const auto find_ancestor = [&]( const Dart& leaf_dart ) {
        return iterateAncestors( leaf_dart, [&]( const Dart& ancestor_global_d ) {
            if( mUnrefinedDarts.at( ancestor_global_d.id() ) )
            {
                global_ancestor_d.emplace( ancestor_global_d );
                return false;
            }
            return true;
        } );
    };

    iterateDartsOfRestrictedCell( *this, leaf_c, dim(), find_ancestor );

    if( not global_ancestor_d ) throw std::runtime_error( "Unable to find an unrefined ancestor dart" );

    return mRanges.toLocalDart( global_ancestor_d.value() );
}

SmallVector<std::variant<bool, size_t>, 3> tpDirectionAlongTPDartPos( const SmallVector<TPCombinatorialMap::TPDartPos, 2>& dart_pos, const util::IndexVec& lengths )
{
    using TPDartPos = TPCombinatorialMap::TPDartPos;
    if( dart_pos.size() == 1 )
    {
        switch( dart_pos.at( 0 ) )
        {
            case TPDartPos::DartPos0:
                return { true, size_t{0} };
            case TPDartPos::DartPos1:
                return { lengths.at( 0 ) - 1, true };
            case TPDartPos::DartPos2:
                return { false, lengths.at( 1 ) - 1 };
            case TPDartPos::DartPos3:
                return { size_t{0}, false };
            default:
                throw std::runtime_error( "Invalid TPDartPos for 2d" );
        }
    }
    else if( dart_pos.size() == 2 )
    {
        const auto dir2d = tpDirectionAlongTPDartPos( { dart_pos.at( 0 ) }, { lengths.at( 0 ), lengths.at( 1 ) } );
        SmallVector<std::variant<bool, size_t>, 3> out;
        switch( dart_pos.at( 1 ) )
        {
            case TPDartPos::DartPos0:
                out = dir2d;
                out.push_back( size_t{0} );
                return out;
            case TPDartPos::DartPos1:
                std::transform( dir2d.begin(), dir2d.end(), std::back_inserter( out ), []( const auto& v ) -> std::variant<bool, size_t> {
                    if( std::holds_alternative<bool>( v ) ) return not std::get<bool>( v );
                    else return std::get<size_t>( v );
                } );
                out.push_back( size_t{0} );
                return out;
            case TPDartPos::DartPos2:
                std::transform( dir2d.begin(),
                                dir2d.end(),
                                lengths.begin(),
                                std::back_inserter( out ),
                                []( const auto& v, const size_t len ) -> std::variant<bool, size_t> {
                                    if( std::holds_alternative<bool>( v ) )
                                        return std::get<bool>( v ) ? size_t{ 0 } : len - 1;
                                    else
                                        return std::get<size_t>( v );
                                } );
                out.push_back( true );
                return out;
            case TPDartPos::DartPos3:
                out = dir2d;
                out.push_back( lengths.at( 2 ) - 1 );
                return out;
            case TPDartPos::DartPos4:
                std::transform( dir2d.begin(),
                                dir2d.end(),
                                lengths.begin(),
                                std::back_inserter( out ),
                                []( const auto& v, const size_t len ) -> std::variant<bool, size_t> {
                                    if( std::holds_alternative<bool>( v ) )
                                        return std::get<bool>( v ) ? len - 1 : size_t{ 0 };
                                    else
                                        return std::get<size_t>( v );
                                } );
                out.push_back( false );
                return out;
            case TPDartPos::DartPos5:
                std::transform( dir2d.begin(), dir2d.end(), std::back_inserter( out ), []( const auto& v ) -> std::variant<bool, size_t> {
                    if( std::holds_alternative<bool>( v ) ) return not std::get<bool>( v );
                    else return std::get<size_t>( v );
                } );
                out.push_back( lengths.at( 2 ) - 1 );
                return out;
            default:
                throw std::runtime_error( "Invalid TPDartPos for 3d" );
        }
    }
 
    throw std::invalid_argument( "No TPDartPos passed to tpDirectionAlongTPDartPos" );
}

bool HierarchicalTPCombinatorialMap::iterateDartLineage( const Dart& global_d,
                                                         const size_t ancestor_or_descendant_level,
                                                         const std::function<bool( const Dart& )>& callback ) const
{
    const auto [dart_level, local_d] = mRanges.toLocalDart( global_d );

    if( dart_level == ancestor_or_descendant_level ) return callback( global_d );
    
    if( dart_level > ancestor_or_descendant_level )
    {
        // Call back on a single ancestor dart
        const size_t darts_per_ancestor_dart =
            std::reduce( std::next( mRefinementRatios.begin(), ancestor_or_descendant_level ),
                         std::next( mRefinementRatios.begin(), dart_level ),
                         1, std::multiplies<>() );

        const FullyUnflattenedDart unflat = unflattenFull( *mRefinementLevels.at( dart_level ), local_d );

        if( checkForNoAncestor( unflat, darts_per_ancestor_dart ) ) return true;

        FullyUnflattenedDart ancestor_unflat( {}, unflat.dart_pos );
        std::transform(
            unflat.unflat_darts.begin(),
            unflat.unflat_darts.end(),
            std::back_inserter( ancestor_unflat.unflat_darts ),
            [&darts_per_ancestor_dart]( const Dart& d ) { return Dart( d.id() / darts_per_ancestor_dart ); } );

        const Dart ancestor_dart = flattenFull( *mRefinementLevels.at( ancestor_or_descendant_level ), ancestor_unflat );
        return callback( mRanges.toGlobalDart( ancestor_or_descendant_level, ancestor_dart ) );
    }
    else
    {
        // Iterate several descendants.
        const size_t darts_per_ancestor_dart =
            std::reduce( std::next( mRefinementRatios.begin(), dart_level ),
                         std::next( mRefinementRatios.begin(), ancestor_or_descendant_level ),
                         1, std::multiplies<>() );

        const FullyUnflattenedDart unflat = unflattenFull( *mRefinementLevels.at( dart_level ), local_d );

        FullyUnflattenedDart start_dart( {}, unflat.dart_pos );
        std::transform(
            unflat.unflat_darts.begin(),
            unflat.unflat_darts.end(),
            std::back_inserter( start_dart.unflat_darts ),
            [&darts_per_ancestor_dart]( const Dart& d ) { return Dart( d.id() * darts_per_ancestor_dart ); } );

        // Add a series of TP indices to the start_dart, flatten, and call back.
        const util::IndexVec lengths( dim(), darts_per_ancestor_dart );
        const SmallVector<std::variant<bool, size_t>, 3> direction = tpDirectionAlongTPDartPos( unflat.dart_pos, lengths );
        const util::IndexVec order = [this]() {
            util::IndexVec order( dim() );
            std::iota( order.begin(), order.end(), 0 );
            return order;
        }();

        bool continue_iter = true;
        util::iterateTensorProduct( lengths, order, direction, [&]( const util::IndexVec& iv ){
            if( not continue_iter ) return;
            FullyUnflattenedDart descendant_unflat = start_dart;
            std::transform( descendant_unflat.unflat_darts.begin(),
                            descendant_unflat.unflat_darts.end(),
                            iv.begin(),
                            descendant_unflat.unflat_darts.begin(),
                            []( const Dart& d, const size_t offset ) { return Dart( d.id() + offset ); } );
            continue_iter = callback( mRanges.toGlobalDart(
                ancestor_or_descendant_level,
                flattenFull( *mRefinementLevels.at( ancestor_or_descendant_level ), descendant_unflat ) ) );
        } );

        return continue_iter;
    }
}

bool HierarchicalTPCombinatorialMap::iterateAncestors( const Dart& global_d,
                                                       const std::function<bool( const Dart& )>& callback ) const
{
    const auto [dart_level, local_d] = mRanges.toLocalDart( global_d );
    for( size_t i = 1; i <= dart_level; i++ )
    {
        const size_t ancestor_level = dart_level - i;
        if( not iterateDartLineage( global_d, ancestor_level, callback ) ) return false;
    }
    return true;
}

bool HierarchicalTPCombinatorialMap::iterateLeafDescendants( const Dart& global_d,
                                                             const std::function<bool( const Dart& )>& callback ) const
{
    const std::function<bool( const Dart& )> recursive_depth_first_iter = [&]( const Dart& global_d ){
        if( mLeafDarts.at( global_d.id() ) )
            return callback( global_d );
        else
        {
            const auto [ level, level_d ] = mRanges.toLocalDart( global_d );
            if( level == mRefinementLevels.size() - 1 )
                throw std::runtime_error( "No descendant leaf dart found" );
            return iterateDartLineage( global_d, level + 1, recursive_depth_first_iter );
        }
    };

    return recursive_depth_first_iter( global_d );
}

bool HierarchicalTPCombinatorialMap::iterateChildren( const Cell& local_cell,
                                                      const size_t cell_level,
                                                      const std::function<bool( const Cell& )>& callback ) const
{
    if ( local_cell.dim() != dim() )
        throw std::runtime_error( "iterateChildren only supports elements." );

    const size_t descendant_level = cell_level + 1;

    if( mRefinementLevels.size() <= descendant_level ) return true;

    // Iterate several descendants.
    const size_t darts_per_ancestor_dart = mRefinementRatios.at( cell_level );

    const FullyUnflattenedDart unflat = unflattenFull( *mRefinementLevels.at( cell_level ), local_cell.dart() );

    FullyUnflattenedDart start_dart( {}, unflat.dart_pos );
    std::transform(
        unflat.unflat_darts.begin(),
        unflat.unflat_darts.end(),
        std::back_inserter( start_dart.unflat_darts ),
        [&darts_per_ancestor_dart]( const Dart& d ) { return Dart( d.id() * darts_per_ancestor_dart ); } );

    // Add a series of TP indices to the start_dart, flatten, and call back.
    const util::IndexVec lengths( dim(), darts_per_ancestor_dart );

    bool continue_iter = true;
    util::iterateTensorProduct( lengths, [&]( const util::IndexVec& iv ){
        if( not continue_iter ) return;
        FullyUnflattenedDart descendant_unflat = start_dart;
        std::transform( descendant_unflat.unflat_darts.begin(),
                        descendant_unflat.unflat_darts.end(),
                        iv.begin(),
                        descendant_unflat.unflat_darts.begin(),
                        []( const Dart& d, const size_t offset ) { return Dart( d.id() + offset ); } );
        continue_iter = callback( topology::Cell( flattenFull( *mRefinementLevels.at( descendant_level ), descendant_unflat ), local_cell.dim() ) );
    } );

    return continue_iter;
}

namespace topology
{
    std::vector<std::vector<Cell>> leafElements( const HierarchicalTPCombinatorialMap& cmap )
    {
        std::vector<std::vector<Cell>> out( cmap.numLevels(), std::vector<Cell>() );
        iterateCellsWhile( cmap, cmap.dim(), [&]( const Cell& c ) {
            const auto [level, dart] = cmap.unrefinedAncestorDartOfCell( c );
            out.at( level ).push_back( Cell( dart, cmap.dim() ) );
            return true;
        } );
        return out;
    }
}