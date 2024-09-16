#include <HierarchicalTPCombinatorialMap.hpp>
#include <GlobalCellMarker.hpp>
#include <IndexOperations.hpp>
#include <ranges>
#include <numeric>

using namespace topology;

HierarchicalTPCombinatorialMap::HierarchicalTPCombinatorialMap(
    const std::vector<std::shared_ptr<const TPCombinatorialMap>>& refinement_levels,
    const std::vector<std::vector<Cell>>& leaf_elements )
    : mRefinementLevels( refinement_levels ),
      mRanges( initializeRanges( refinement_levels ) ),
      mLeafDarts( mRanges.maxDartId(), false ),
      mUnrefinedDarts( mRanges.maxDartId(), false )
{
    // TODO: For now I am assuming dyadic refinement only.  This checks that is the case, though imperfectly.
    for( size_t level = 1; level < refinement_levels.size(); level++ )
    {
        if( topology::cellCount( *refinement_levels.at( level - 1 ), dim() ) * 4 !=
            topology::cellCount( *refinement_levels.at( level ), dim() ) )
            throw std::runtime_error( "Hierarchical only supports dyadic refinement currently.  Input refinement levels do not form a dyadic refinement." );
    }

    // Mark unrefined and leaf darts
    for( size_t level = 0; level < refinement_levels.size(); level++ )
    {
        for( const auto& leaf_elem : leaf_elements.at( level ) )
        {
            iterateDartsOfCell( *refinement_levels.at( level ), leaf_elem, [&]( const Dart& level_d ) {
                mUnrefinedDarts.at( mRanges.toGlobalDart( level, level_d ).id() ) = true;
                for( size_t cell_dim = 1; cell_dim < dim(); cell_dim++ )
                {
                    const Cell c( level_d, cell_dim );
                    iterateDartsOfCell( *refinement_levels.at( level ), c, [&]( const Dart& adj_level_d ) {
                        const Dart adj_global_d = mRanges.toGlobalDart( level, adj_level_d );
                        mLeafDarts.at( adj_global_d.id() ) = true;

                        // Unmark ancestors that are marked as leaf.
                        iterateAncestors( adj_global_d, [&]( const Dart& ancestor_global_d ) {
                            mLeafDarts.at( ancestor_global_d.id() ) = false;
                            return true;
                        } );
                        return true;
                    } );
                }
                return true;
            } );
        }
    }

    // Build phi1 and phi-1 ops
    for( size_t level = 0; level < refinement_levels.size(); level++ )
    {
        const TPCombinatorialMap& level_cmap = *refinement_levels.at( level );
        for( const auto& leaf_elem : leaf_elements.at( level ) )
        {
            const bool all_darts_are_leaves = iterateDartsOfCell( level_cmap, leaf_elem, [&]( const Dart& d ) {
                const Dart global_d = mRanges.toGlobalDart( level, d );
                return mLeafDarts.at( global_d.id() );
            } );

            if( not all_darts_are_leaves )
            {
                std::optional<Dart> first_leaf;
                std::optional<Dart> previous_leaf;
                const Dart global_leaf_elem_dart = mRanges.toGlobalDart( level, leaf_elem.dart() );
                Dart d = mRanges.toGlobalDart( level, leaf_elem.dart() );
                do
                {
                    iterateLeafDescendants( d, [&]( const Dart& descendant ) {
                        if( previous_leaf )
                        {
                            mPhiOnes.insert( { previous_leaf.value(), descendant } );
                            mPhiMinusOnes.insert( { descendant, previous_leaf.value() } );
                        }
                        previous_leaf.emplace( descendant );
                        if( not first_leaf ) first_leaf.emplace( descendant );
                        return true;
                    } );
                    d = mRanges.toGlobalDart( level, topology::phi( level_cmap, 1, mRanges.toLocalDart( d ).second ).value() );
                } while ( d != global_leaf_elem_dart );
                
                if( first_leaf )
                {
                    mPhiOnes.insert( { previous_leaf.value(), first_leaf.value() } );
                    mPhiMinusOnes.insert( { first_leaf.value(), previous_leaf.value() } );
                }
            }
        }
    }
}

std::optional<Dart> HierarchicalTPCombinatorialMap::phi( const int i, const Dart& d ) const
{
    const auto [level, level_d ] = mRanges.toLocalDart( d );
    return topology::phi( *mRefinementLevels.at( level ), i, level_d ).and_then( [&]( const Dart& level_phi ) -> std::optional<Dart> {
        const Dart level_global_phi = mRanges.toGlobalDart( level, level_phi );

        if( i > 1 or mLeafDarts.at( level_global_phi.id() ) ) return level_global_phi;

        const std::map<Dart, Dart>& phi_map = i == 1 ? mPhiOnes : mPhiMinusOnes;
        const auto it = phi_map.find( d );
        if( it != phi_map.end() ) return it->second;
        else throw std::runtime_error( "Phi 1 or -1 operation unfound!" );
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

std::pair<size_t, Dart> HierarchicalTPCombinatorialMap::unrefinedAncestorDart( const Dart& leaf_d ) const
{
    if( not mLeafDarts.at( leaf_d.id() ) )
        throw std::runtime_error( "unrefinedAncestorDart requires a leaf dart as input" );

    if( mUnrefinedDarts.at( leaf_d.id() ) ) return mRanges.toLocalDart( leaf_d );

    std::optional<Dart> global_ancestor_d;

    iterateAncestors( leaf_d, [&]( const Dart& ancestor_global_d ) {
        if( mUnrefinedDarts.at( ancestor_global_d.id() ) )
        {
            global_ancestor_d.emplace( ancestor_global_d );
            return false;
        }
        return true;
    } );

    if( not global_ancestor_d ) throw std::runtime_error( "Unable to find an unrefined ancestor dart" );

    return mRanges.toLocalDart( global_ancestor_d.value() );
}

SmallVector<std::variant<bool, size_t>, 3> tpDirectionAlongTPDartPos( const SmallVector<TPCombinatorialMap::TPDartPos, 2>& dart_pos, const util::IndexVec& lengths )
{
    using TPDartPos = TPCombinatorialMap::TPDartPos;
    SmallVector<std::variant<bool, size_t>, 3> out;
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
    throw std::runtime_error( "Hierarchical not implemented yet for 3d" );
}

bool checkForNoAncestor( const FullyUnflattenedDart& unflat, const size_t n_darts_per_ancestor )
{
    if( unflat.unflat_darts.size() == 3 ) throw std::runtime_error( "Hierarchical not implemented yet for 3d" );

    // 2D case

    switch( unflat.dart_pos.at( 0 ) )
    {
        case TPCombinatorialMap::TPDartPos::DartPos0:
            return unflat.unflat_darts.at( 1 ).id() % n_darts_per_ancestor != 0;
        case TPCombinatorialMap::TPDartPos::DartPos1:
            return ( unflat.unflat_darts.at( 0 ).id() + 1 ) % n_darts_per_ancestor != 0;
        case TPCombinatorialMap::TPDartPos::DartPos2:
            return ( unflat.unflat_darts.at( 1 ).id() + 1 ) % n_darts_per_ancestor != 0;
        case TPCombinatorialMap::TPDartPos::DartPos3:
            return unflat.unflat_darts.at( 0 ).id() % n_darts_per_ancestor != 0;
        default:
            throw std::runtime_error( "Invalid TPDartPos for 2d" );
    }
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
        // TODO: Assumes dyadic refinement, which is currently checked for in constructor.
        const size_t darts_per_ancestor_dart = 2 * ( dart_level - ancestor_or_descendant_level );

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
        // TODO: Assumes dyadic refinement, which is currently checked for in constructor.
        const size_t darts_per_ancestor_dart = 2 * ( ancestor_or_descendant_level - dart_level );

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

    // Iterate several descendants.
    // TODO: Assumes dyadic refinement, which is currently checked for in constructor.
    constexpr size_t darts_per_ancestor_dart = 2;

    if( mRefinementLevels.size() <= descendant_level ) return true;

    const FullyUnflattenedDart unflat = unflattenFull( *mRefinementLevels.at( cell_level ), local_cell.dart() );

    FullyUnflattenedDart start_dart( {}, unflat.dart_pos );
    std::transform(
        unflat.unflat_darts.begin(),
        unflat.unflat_darts.end(),
        std::back_inserter( start_dart.unflat_darts ),
        []( const Dart& d ) { return Dart( d.id() * darts_per_ancestor_dart ); } );

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