#include <HierarchicalMultiPatchCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>
#include <numeric>

using namespace topology;

std::vector<std::shared_ptr<const HierarchicalTPCombinatorialMap>>
    initializeConstituents( const std::vector<std::shared_ptr<const MultiPatchCombinatorialMap>>& refinement_levels,
                            const std::vector<std::vector<Cell>>& leaf_elements )
{
    if( refinement_levels.empty() )
        throw std::invalid_argument( "HierarchicalMultiPatchCombinatorialMap requires at least one refinement level." );
    if( leaf_elements.size() != refinement_levels.size() )
        throw std::invalid_argument( "Incorrect number of leaf element levels given to HierarchicalMultiPatchCombinatorialMap." );

    const uint dim = refinement_levels.front()->dim();

    //////////////////////////////////////////////////////////////////////////////
    /// initialize mutable constituents

    std::vector<std::vector<std::shared_ptr<const TPCombinatorialMap>>> constituent_levels(
        refinement_levels.front()->constituents().size() );
    for( auto& level : constituent_levels )
        level.reserve( refinement_levels.size() );

    std::vector<std::vector<std::vector<Cell>>> constituent_leaf_elements(
        refinement_levels.front()->constituents().size(), std::vector<std::vector<Cell>>( refinement_levels.size() ) );

    for( size_t level_ii = 0; level_ii < refinement_levels.size(); level_ii++ )
    {
        const auto& level = refinement_levels.at( level_ii );
        if( level->constituents().size() != constituent_levels.size() )
            throw std::invalid_argument( "All refinement levels must have the same number of constituents." );

        for( size_t patch_ii = 0; patch_ii < level->constituents().size(); patch_ii++ )
        {
            constituent_levels.at( patch_ii ).push_back( level->constituents().at( patch_ii ) );
        }

        for( const auto& elem : leaf_elements.at( level_ii ) )
        {
            const auto [patch_ii, local_d] = refinement_levels.at( level_ii )->toLocalDart( elem.dart() );
            constituent_leaf_elements.at( patch_ii ).at( level_ii ).push_back( Cell( local_d, elem.dim() ) );
        }
    }

    std::vector<std::unique_ptr<MutableHierarchicalTPCombinatorialMap>> mutable_constituents;

    for( size_t patch_ii = 0; patch_ii < constituent_levels.size(); patch_ii++ )
    {
        mutable_constituents.push_back( std::make_unique<MutableHierarchicalTPCombinatorialMap>(
            constituent_levels.at( patch_ii ), constituent_leaf_elements.at( patch_ii ) ) );
    }

    /////////////////////////////////////////////////////////////////////////////////
    /// Find the darts that need to be refined because of the multi-patch connections

    std::set<std::pair<size_t, Dart>> not_leaf;

    const auto mark_as_leaf = [&]( const size_t patch_ii, const size_t level_ii, const Dart& constituent_level_d ) {
        auto& mutable_constituent = *mutable_constituents.at( patch_ii );
        const Dart constituent_dart = mutable_constituent.dartRanges().toGlobalDart( level_ii, constituent_level_d );
        mutable_constituent.setLeaf( constituent_dart, true );

        // Unmark ancestors that are marked as leaf.
        mutable_constituent.iterateAncestors( constituent_dart, [&]( const Dart& ancestor_global_d ) {
            not_leaf.emplace( patch_ii, ancestor_global_d );
            return true;
        } );
    };

    const auto has_leaf_ancestor = [&]( const size_t patch_ii, const size_t level_ii, const Dart& constituent_level_d ) {
        const auto& mutable_constituent = *mutable_constituents.at( patch_ii );
        const Dart constituent_dart = mutable_constituent.dartRanges().toGlobalDart( level_ii, constituent_level_d );
        return not mutable_constituent.iterateAncestors( constituent_dart, [&]( const Dart& ancestor_global_d ) {
            return not mutable_constituent.isLeaf( ancestor_global_d.id() );
        } );
    };

    // Mark unrefined and leaf darts between boundaries.  Don't worry about level 0, those will be handled by the constituents.
    for( size_t level = 1; level < refinement_levels.size(); level++ )
    {
        for( const auto& leaf_elem : leaf_elements.at( level ) )
        {
            iterateDartsOfCell( *refinement_levels.at( level ), leaf_elem, [&]( const Dart& level_d ) {
                // We only need to do something here if this is a boundary dart in the constituent, but not in the multipatch.
                const auto [patch_ii, constituent_level_d] = refinement_levels.at( level )->toLocalDart( level_d );
                const auto maybe_constituent_level_phi = phi( *refinement_levels.at( level )->constituents().at( patch_ii ), dim, constituent_level_d );
                if( maybe_constituent_level_phi.has_value() ) return true; // This is not a boundary dart in the constituent.
                const auto maybe_phi = phi( *refinement_levels.at( level ), dim, level_d );
                if( not maybe_phi.has_value() ) return true; // This is a boundary dart in the multipatch.

                mark_as_leaf( patch_ii, level, constituent_level_d );
                const auto [other_patch_ii, other_constituent_level_d] = refinement_levels.at( level )->toLocalDart( maybe_phi.value() );
                mark_as_leaf( other_patch_ii, level, other_constituent_level_d );
                if( dim == 3 )
                {
                    iterateDartsOfCell( *refinement_levels.at( level ), Edge( level_d ), [&]( const Dart& other_level_d ) {
                        const auto [other_patch_ii, other_constituent_level_d] = refinement_levels.at( level )->toLocalDart( other_level_d );
                        if( has_leaf_ancestor( other_patch_ii, level, other_constituent_level_d ) )
                        {
                            mark_as_leaf( other_patch_ii, level, other_constituent_level_d );
                            const auto maybe_other_phi = phi( *refinement_levels.at( level ), 3, other_level_d );
                            if( maybe_other_phi.has_value() )
                            {
                                const auto [other_patch_ii, other_constituent_level_d] = refinement_levels.at( level )->toLocalDart( maybe_other_phi.value() );
                                mark_as_leaf( other_patch_ii, level, other_constituent_level_d );
                            }
                        }
                        return true;
                    } );
                }

                return true;
            } );
        }
    }

    for( const std::pair<size_t, Dart>& d : not_leaf )
    {
        auto& mutable_constituent = *mutable_constituents.at( d.first );
        mutable_constituent.setLeaf( d.second, false );
    }

    // Build phi1 and phi-1 ops
    for( size_t level = 0; level < refinement_levels.size() - 1; level++ )
    {
        for( const auto& leaf_elem : leaf_elements.at( level ) )
        {
            const MultiPatchCombinatorialMap& refinement_level = *refinement_levels.at( level );
            const auto [patch_ii, patch_level_dart] = refinement_level.toLocalDart( leaf_elem.dart() );

            const topology::Cell patch_leaf_elem( patch_level_dart, leaf_elem.dim() );

            const auto phi_chain = [&level,&refinement_level,&mutable_constituents,&dim]( const size_t patch_ii, const topology::Face& leaf_face ){
                const TPCombinatorialMap& level_patch_cmap = *refinement_level.constituents().at( patch_ii );
                auto& mutable_constituent = *mutable_constituents.at( patch_ii );

                const bool all_darts_are_leaves = iterateDartsOfRestrictedCell( level_patch_cmap, leaf_face, 3, [&]( const Dart& d ) {
                    const Dart global_d = mutable_constituent.dartRanges().toGlobalDart( level, d );
                    return mutable_constituent.isLeaf( global_d.id() );
                } );

                const bool all_darts_are_not_leaves = iterateDartsOfRestrictedCell( level_patch_cmap, leaf_face, 3, [&]( const Dart& d ) {
                    const Dart global_d = mutable_constituent.dartRanges().toGlobalDart( level, d );
                    return not mutable_constituent.isLeaf( global_d.id() );
                } );

                if( not all_darts_are_leaves and not ( dim == 3 and all_darts_are_not_leaves ) )
                {
                    std::optional<Dart> first_leaf;
                    std::optional<Dart> previous_leaf;
                    const Dart global_leaf_face_dart = mutable_constituent.dartRanges().toGlobalDart( level, leaf_face.dart() );
                    Dart d = global_leaf_face_dart;
                    do
                    {
                        mutable_constituent.iterateLeafDescendants( d, [&]( const Dart& descendant ) {
                            if( previous_leaf )
                            {
                                mutable_constituent.setPhi( previous_leaf.value(), descendant );
                            }
                            previous_leaf.emplace( descendant );
                            if( not first_leaf ) first_leaf.emplace( descendant );
                            return true;
                        } );
                        d = mutable_constituent.dartRanges().toGlobalDart( level, topology::phi( level_patch_cmap, 1, mutable_constituent.dartRanges().toLocalDart( d ).second ).value() );
                    } while ( d != global_leaf_face_dart );
 
                    if( first_leaf )
                    {
                        mutable_constituent.setPhi( previous_leaf.value(), first_leaf.value() );
                    }
                }
            };

            // In 3d we need to do this for the faces of the leaf elements, in 2d, the elements themselves.
            const TPCombinatorialMap& level_patch_cmap = *refinement_level.constituents().at( patch_ii );
            iterateAdjacentCells( level_patch_cmap, patch_leaf_elem, 2, [&]( const topology::Face& patch_leaf_face ) {
                phi_chain( patch_ii, patch_leaf_face );
                if( dim == 3 )
                {
                    const Dart multipatch_leaf_dart = refinement_level.toGlobalDart( patch_ii, patch_leaf_face.dart() );
                    const auto maybe_phi3 = phi( refinement_level, 3, multipatch_leaf_dart );
                    if( maybe_phi3.has_value() )
                    {
                        const auto [other_patch_ii, other_patch_level_dart] = refinement_level.toLocalDart( maybe_phi3.value() );
                        phi_chain( other_patch_ii, other_patch_level_dart );
                    }
                }
                return true;
            } );
        }
    }

    /////////////////////////////////////////////////////////////////////////////////
    /// Make constituents immutable
    std::vector<std::shared_ptr<const HierarchicalTPCombinatorialMap>> constituents;
    for( const auto& con : mutable_constituents )
    {
        constituents.push_back( std::make_shared<const HierarchicalTPCombinatorialMap>( con->asImmutable() ) );
    }
    return constituents;
}

HierarchicalMultiPatchCombinatorialMap::HierarchicalMultiPatchCombinatorialMap(
    const std::vector<std::shared_ptr<const MultiPatchCombinatorialMap>>& refinement_levels,
    const std::vector<std::vector<Cell>>& leaf_elements ) :
    mRefinementLevels( refinement_levels ), mConstituents( initializeConstituents( refinement_levels, leaf_elements ) ), mRanges( initializeRanges( constituents() ) )
{}

std::optional<Dart> HierarchicalMultiPatchCombinatorialMap::phi( const int i, const Dart& d ) const
{
    // Convert to a constituent dart, then do the phi there.
    // If the phi is a nothing, and i == dim, then check the phi(dim) in the multipatch of the correct level.

    const auto [patch_ii, local_d] = mRanges.toLocalDart( d );
    return mConstituents.at( patch_ii )->phi( i, local_d ).and_then( [&]( const Dart& local_phi ) -> std::optional<Dart> {
        return mRanges.toGlobalDart( patch_ii, local_phi );
    } ).or_else( [&]() -> std::optional<Dart> {
        if( i == (int)dim() )
        {
            const auto [level_ii, level_constituent_d] = mConstituents.at( patch_ii )->dartRanges().toLocalDart( local_d );
            const Dart level_d = mRefinementLevels.at( level_ii )->toGlobalDart( patch_ii, level_constituent_d );
            return mRefinementLevels.at( level_ii )->phi( i, level_d ).and_then( [&]( const Dart& level_phi ) -> std::optional<Dart> {
                const auto [patch_ii, level_constituent_phi] = mRefinementLevels.at( level_ii )->toLocalDart( level_phi );
                const Dart local_phi = mConstituents.at( patch_ii )->dartRanges().toGlobalDart( level_ii, level_constituent_phi );
                return mRanges.toGlobalDart( patch_ii, local_phi );
            } );
        }
        else return std::nullopt;
    } );
}

Dart::IndexType HierarchicalMultiPatchCombinatorialMap::maxDartId() const
{
    return mRanges.maxDartId();
}

uint HierarchicalMultiPatchCombinatorialMap::dim() const
{
    return mRefinementLevels.front()->dim();
}

bool HierarchicalMultiPatchCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    for( size_t patch_ii = 0; patch_ii < mConstituents.size(); patch_ii++ )
    {
        const auto& constituent = *mConstituents.at( patch_ii );
        const bool continue_iter = topology::iterateDartsWhile( constituent, [&]( const Dart& constituent_d ) {
            return callback( mRanges.toGlobalDart( patch_ii, constituent_d ) );
        } );
        if( not continue_iter ) return false;
    }
    return true;
}

bool HierarchicalMultiPatchCombinatorialMap::iterateCellsWhile( const uint cell_dim,
                            const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim == dim() )
    {
        for( size_t patch_ii = 0; patch_ii < mConstituents.size(); patch_ii++ )
        {
            const auto& constituent = *mConstituents.at( patch_ii );
            const bool continue_iter = topology::iterateCellsWhile( constituent, cell_dim, [&]( const Cell& constituent_c ) {
                return callback( Cell( mRanges.toGlobalDart( patch_ii, constituent_c.dart() ), constituent_c.dim() ) );
            } );
            if( not continue_iter ) return false;
        }
        return true;
    }
    else
    {
        GlobalCellMarker m( *this, cell_dim );
        return iterateDartsWhile( [&]( const Dart& d ) {
            const Cell c( d, cell_dim );
            if( not m.isMarked( c ) )
            {
                m.mark( *this, c );
                return callback( c );
            }
            return true;
        } );
    }
}

std::optional<IndexingFunc> HierarchicalMultiPatchCombinatorialMap::indexing( const uint ) const
{
    return std::nullopt;
}

std::optional<size_t> HierarchicalMultiPatchCombinatorialMap::cellCount( const uint cell_dim ) const
{
    if( cell_dim == dim() )
        return std::transform_reduce( mConstituents.begin(),
                                      mConstituents.end(),
                                      0,
                                      std::plus<>(),
                                      [&cell_dim]( const auto& m ) {
                                          return topology::cellCount( *m, cell_dim );
                                      } );
    else return std::nullopt;
}

std::pair<size_t, Dart> HierarchicalMultiPatchCombinatorialMap::unrefinedAncestorDart( const Dart& leaf_d ) const
{
    const auto [patch_ii, local_d] = mRanges.toLocalDart( leaf_d );
    const auto [level, level_d] = mConstituents.at( patch_ii )->unrefinedAncestorDart( local_d );
    return { level, mRefinementLevels.at( level )->toGlobalDart( patch_ii, level_d ) };
}

bool HierarchicalMultiPatchCombinatorialMap::iterateChildren( const Cell& local_cell,
                            const size_t cell_level,
                            const std::function<bool( const Cell& )>& callback ) const
{
    if( local_cell.dim() != dim() )
        throw std::runtime_error( "iterateChildren only supports elements." );

    const auto [patch_ii, local_d] = mRefinementLevels.at( cell_level )->toLocalDart( local_cell.dart() );
    return mConstituents.at( patch_ii )->iterateChildren( Cell( local_d, local_cell.dim() ), cell_level, [&]( const Cell& child_constituent_c ) {
        return callback( Cell( mRefinementLevels.at( cell_level + 1 )->toGlobalDart( patch_ii, child_constituent_c.dart() ), child_constituent_c.dim() ) );
    } );
}

namespace topology
{
    std::vector<std::vector<Cell>> leafElements( const HierarchicalMultiPatchCombinatorialMap& cmap )
    {
        std::vector<std::vector<Cell>> leaf_elements( cmap.numLevels() );
        for( size_t patch_ii = 0; patch_ii < cmap.constituents().size(); patch_ii++ )
        {
            const auto& constituent = cmap.constituents().at( patch_ii );
            const std::vector<std::vector<Cell>> constituent_leaf_elements = leafElements( *constituent );
            for( size_t level = 0; level < constituent_leaf_elements.size(); level++ )
            {
                const auto& level_elems = constituent_leaf_elements.at( level );
                std::transform(
                    level_elems.begin(),
                    level_elems.end(),
                    std::back_inserter( leaf_elements.at( level ) ),
                    [&]( const Cell& c ) {
                        return Cell( cmap.refinementLevels().at( level )->toGlobalDart( patch_ii, c.dart() ), c.dim() );
                    } );
            }
        }
        return leaf_elements;
    }
}; // namespace topology