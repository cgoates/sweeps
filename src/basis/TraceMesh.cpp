#include <TraceMesh.hpp>
#include <BasisComplex.hpp>
#include <CombinatorialMapMethods.hpp>
#include <HierarchicalMultiPatchCombinatorialMap.hpp>
#include <IndexOperations.hpp>
#include <ParametricAtlas.hpp>
#include <TPCombinatorialMap.hpp>
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <set>
#include <stdexcept>

namespace basis
{
    namespace
    {
        using TPPermutation = topology::MultiPatchCombinatorialMap::TPPermutation;
        using SideIndex = std::vector<size_t>;

        struct LevelTraceMesh
        {
            size_t level = 0;
            const topology::MultiPatchCombinatorialMap* cmap = nullptr;
            std::vector<topology::Cell> leaf_elements;
            std::function<topology::Cell( const topology::Cell& )> to_space_cell;
        };

        void checkTopDimensionalSpace( const SplineSpace& ss )
        {
            const uint dim = ss.basisComplex().parametricAtlas().cmap().dim();
            if( dim < 2 or dim > 3 )
                throw std::runtime_error( "Trace mesh enumeration is implemented for 2D and 3D spaces." );
        }

        std::vector<LevelTraceMesh> levelTraceMeshes( const SplineSpace& ss )
        {
            checkTopDimensionalSpace( ss );

            const auto* mp_cmap =
                dynamic_cast<const topology::MultiPatchCombinatorialMap*>(
                    &ss.basisComplex().parametricAtlas().cmap() );
            if( mp_cmap != nullptr )
            {
                std::vector<topology::Cell> cells;
                topology::iterateCellsWhile( *mp_cmap, mp_cmap->dim(), [&]( const topology::Cell& c ) {
                    cells.push_back( c );
                    return true;
                } );

                return { LevelTraceMesh{ 0, mp_cmap, cells, []( const topology::Cell& c ) { return c; } } };
            }

            const auto* hmp_cmap =
                dynamic_cast<const topology::HierarchicalMultiPatchCombinatorialMap*>(
                    &ss.basisComplex().parametricAtlas().cmap() );
            if( hmp_cmap == nullptr )
                throw std::runtime_error( "Trace mesh enumeration requires a multipatch or hierarchical multipatch space." );

            const std::vector<std::vector<topology::Cell>> leaves = topology::leafElements( *hmp_cmap );
            std::vector<LevelTraceMesh> out;
            out.reserve( hmp_cmap->numLevels() );
            for( size_t level = 0; level < hmp_cmap->numLevels(); level++ )
            {
                out.push_back( LevelTraceMesh{
                    level,
                    hmp_cmap->refinementLevels().at( level ).get(),
                    leaves.at( level ),
                    [hmp_cmap, level]( const topology::Cell& c ) {
                        return topology::Cell( hmp_cmap->toGlobalDart( level, c.dart() ), c.dim() );
                    } } );
            }
            return out;
        }

        std::vector<size_t> sideLengths( const topology::TPCombinatorialMap& cmap, const ElementSide& side )
        {
            const SmallVector<std::shared_ptr<const topology::CombinatorialMap1d>, 3> components =
                topology::tensorProductComponentCMaps( cmap );

            std::vector<size_t> out;
            for( size_t axis = 0; axis < cmap.dim(); axis++ )
            {
                if( axis == side.axis ) continue;
                const auto& component = components.at( axis );
                out.push_back( topology::cellCount( *component, component->dim() ) );
            }
            return out;
        }

        size_t normalLength( const topology::TPCombinatorialMap& cmap, const ElementSide& side )
        {
            const SmallVector<std::shared_ptr<const topology::CombinatorialMap1d>, 3> components =
                topology::tensorProductComponentCMaps( cmap );
            const auto& component = components.at( side.axis );
            return topology::cellCount( *component, component->dim() );
        }

        topology::Cell localCellOnSide( const topology::TPCombinatorialMap& cmap,
                                        const ElementSide& side,
                                        const SideIndex& side_index )
        {
            topology::FullyUnflattenedDart unflat;
            const size_t n_length = normalLength( cmap, side );
            for( size_t axis = 0, side_axis = 0; axis < cmap.dim(); axis++ )
            {
                if( axis == side.axis )
                    unflat.unflat_darts.push_back( topology::Dart( side.lower ? 0 : n_length - 1 ) );
                else
                    unflat.unflat_darts.push_back( topology::Dart( side_index.at( side_axis++ ) ) );
            }
            unflat.dart_pos =
                SmallVector<topology::TPCombinatorialMap::TPDartPos, 2>(
                    cmap.dim() - 1, topology::TPCombinatorialMap::TPDartPos::DartPos0 );
            return topology::Cell( topology::flattenFull( cmap, unflat ), cmap.dim() );
        }

        topology::Cell levelCellOnSide( const topology::MultiPatchCombinatorialMap& cmap,
                                        const size_t patch_id,
                                        const ElementSide& side,
                                        const SideIndex& side_index )
        {
            const topology::Cell local_cell = localCellOnSide( *cmap.constituents().at( patch_id ), side, side_index );
            return topology::Cell( cmap.toGlobalDart( patch_id, local_cell.dart() ), cmap.dim() );
        }

        SideIndex sideIndex( const topology::TPCombinatorialMap& cmap,
                             const topology::Cell& local_cell,
                             const ElementSide& side )
        {
            const topology::FullyUnflattenedDart unflat = topology::unflattenFull( cmap, local_cell.dart() );
            SideIndex out;
            for( size_t axis = 0; axis < cmap.dim(); axis++ )
            {
                if( axis != side.axis ) out.push_back( unflat.unflat_darts.at( axis ).id() );
            }
            return out;
        }

        std::set<SideIndex> sideLeafIndices( const topology::MultiPatchCombinatorialMap& cmap,
                                             const std::vector<topology::Cell>& leaf_elements,
                                             const size_t patch_id,
                                             const ElementSide& side )
        {
            const auto& patch_cmap = *cmap.constituents().at( patch_id );
            const size_t n_length = normalLength( patch_cmap, side );

            std::set<SideIndex> out;
            for( const topology::Cell& leaf : leaf_elements )
            {
                if( leaf.dim() != cmap.dim() ) continue;
                const auto [leaf_patch, local_d] = cmap.toLocalDart( leaf.dart() );
                if( leaf_patch != patch_id ) continue;

                const topology::FullyUnflattenedDart unflat = topology::unflattenFull( patch_cmap, local_d );
                const size_t normal_index = unflat.unflat_darts.at( side.axis ).id();
                const bool on_side = side.lower ? normal_index == 0 : normal_index + 1 == n_length;
                if( not on_side ) continue;

                out.insert( sideIndex( patch_cmap, topology::Cell( local_d, leaf.dim() ), side ) );
            }
            return out;
        }

        std::set<std::pair<size_t, topology::Cell>>
            patchLeafSet( const topology::MultiPatchCombinatorialMap& cmap,
                          const std::vector<topology::Cell>& leaf_elements )
        {
            std::set<std::pair<size_t, topology::Cell>> out;
            for( const topology::Cell& leaf : leaf_elements )
            {
                const auto [patch_id, local_d] = cmap.toLocalDart( leaf.dart() );
                out.emplace( patch_id, topology::Cell( local_d, leaf.dim() ) );
            }
            return out;
        }

        bool traceHasBrokenRows( const SplineSpace& ss,
                                 const TraceMeshSide& first,
                                 const TraceMeshSide& second,
                                 const double row_tol )
        {
            const TraceSideData first_trace = traceSideData( ss, first.element, first.side, row_tol );
            const TraceSideData second_trace = traceSideData( ss, second.element, second.side, row_tol );

            std::set<FunctionId> first_ids( first_trace.connectivity.begin(), first_trace.connectivity.end() );
            std::set<FunctionId> second_ids( second_trace.connectivity.begin(), second_trace.connectivity.end() );
            std::vector<FunctionId> shared;
            std::set_intersection( first_ids.begin(),
                                   first_ids.end(),
                                   second_ids.begin(),
                                   second_ids.end(),
                                   std::back_inserter( shared ) );
            return shared.size() < std::min( first_ids.size(), second_ids.size() );
        }
    }

    std::vector<size_t> permuteTraceSideIndex( const std::vector<size_t>& index,
                                               const std::vector<size_t>& lengths,
                                               const TPPermutation permutation )
    {
        const auto flip = [&]( const size_t axis ) { return lengths.at( axis ) - index.at( axis ) - 1; };

        if( permutation == TPPermutation::Flip1d )
        {
            if( index.size() != 1 )
                throw std::runtime_error( "Flip1d trace permutation requires a one-dimensional side index." );
            return { flip( 0 ) };
        }

        if( index.size() != 2 )
            throw std::runtime_error( "3D trace permutation requires a two-dimensional side index." );

        switch( permutation )
        {
            case TPPermutation::ZeroToZero: return { flip( 0 ), index.at( 1 ) };
            case TPPermutation::ZeroToOne: return { flip( 1 ), flip( 0 ) };
            case TPPermutation::ZeroToTwo: return { index.at( 0 ), flip( 1 ) };
            case TPPermutation::ZeroToThree: return { index.at( 1 ), index.at( 0 ) };
            case TPPermutation::Flip1d: break;
        }
        throw std::runtime_error( "Unknown trace permutation." );
    }

    Eigen::VectorXd permuteTraceSidePoint( const Eigen::VectorXd& point, const TPPermutation permutation )
    {
        if( permutation == TPPermutation::Flip1d )
        {
            if( point.size() != 1 )
                throw std::runtime_error( "Flip1d trace permutation requires a one-dimensional side point." );
            return Eigen::Vector<double, 1>( 1.0 - point( 0 ) );
        }

        if( point.size() != 2 )
            throw std::runtime_error( "3D trace permutation requires a two-dimensional side point." );

        Eigen::Vector2d out;
        switch( permutation )
        {
            case TPPermutation::ZeroToZero: out << 1.0 - point( 0 ), point( 1 ); return out;
            case TPPermutation::ZeroToOne: out << 1.0 - point( 1 ), 1.0 - point( 0 ); return out;
            case TPPermutation::ZeroToTwo: out << point( 0 ), 1.0 - point( 1 ); return out;
            case TPPermutation::ZeroToThree: out << point( 1 ), point( 0 ); return out;
            case TPPermutation::Flip1d: break;
        }
        throw std::runtime_error( "Unknown trace permutation." );
    }

    std::vector<TraceMeshInterface> boundaryTraceMeshInterfaces( const SplineSpace& ss )
    {
        std::vector<TraceMeshInterface> out;
        for( const LevelTraceMesh& level : levelTraceMeshes( ss ) )
        {
            for( size_t patch_id = 0; patch_id < level.cmap->constituents().size(); patch_id++ )
            {
                for( size_t side_id = 0; side_id < 2 * level.cmap->dim(); side_id++ )
                {
                    if( level.cmap->connections().contains( { patch_id, side_id } ) ) continue;

                    const ElementSide side = elementSideFromId( side_id );
                    for( const SideIndex& index : sideLeafIndices( *level.cmap, level.leaf_elements, patch_id, side ) )
                    {
                        const topology::Cell level_cell = levelCellOnSide( *level.cmap, patch_id, side, index );
                        out.push_back( { TraceMeshEntityType::Boundary,
                                         TraceMeshSide{ level.to_space_cell( level_cell ), side, patch_id, level.level },
                                         std::nullopt,
                                         std::nullopt } );
                    }
                }
            }
        }
        return out;
    }

    std::vector<TraceMeshInterface> patchTraceMeshInterfaces( const SplineSpace& ss )
    {
        std::vector<TraceMeshInterface> out;
        for( const LevelTraceMesh& level : levelTraceMeshes( ss ) )
        {
            for( const auto& [first, connection] : level.cmap->connections() )
            {
                const topology::MultiPatchCombinatorialMap::ConstituentSide& second = connection.second;
                if( not( first < second ) ) continue;

                const ElementSide first_side = elementSideFromId( first.side_id );
                const ElementSide second_side = elementSideFromId( second.side_id );
                const auto& first_cmap = *level.cmap->constituents().at( first.constituent_id );
                const std::vector<size_t> first_lengths = sideLengths( first_cmap, first_side );

                for( const SideIndex& first_index :
                     sideLeafIndices( *level.cmap, level.leaf_elements, first.constituent_id, first_side ) )
                {
                    const SideIndex second_index =
                        permuteTraceSideIndex( first_index, first_lengths, connection.first );
                    const topology::Cell first_cell =
                        levelCellOnSide( *level.cmap, first.constituent_id, first_side, first_index );
                    const topology::Cell second_cell =
                        levelCellOnSide( *level.cmap, second.constituent_id, second_side, second_index );

                    out.push_back(
                        { TraceMeshEntityType::PatchInterface,
                          TraceMeshSide{ level.to_space_cell( first_cell ),
                                         first_side,
                                         first.constituent_id,
                                         level.level },
                          TraceMeshSide{ level.to_space_cell( second_cell ),
                                         second_side,
                                         second.constituent_id,
                                         level.level },
                          connection.first } );
                }
            }
        }
        return out;
    }

    std::vector<TraceMeshInterface> interiorTraceMeshInterfaces( const SplineSpace& ss,
                                                                 const InteriorTraceMode mode,
                                                                 const double row_tol )
    {
        std::vector<TraceMeshInterface> out;
        for( const LevelTraceMesh& level : levelTraceMeshes( ss ) )
        {
            const auto leaves = patchLeafSet( *level.cmap, level.leaf_elements );
            std::set<std::pair<topology::Cell, topology::Cell>> seen_pairs;

            for( const topology::Cell& level_leaf : level.leaf_elements )
            {
                const auto [patch_id, local_d] = level.cmap->toLocalDart( level_leaf.dart() );
                const auto& patch_cmap = *level.cmap->constituents().at( patch_id );
                const topology::Cell local_leaf( local_d, level_leaf.dim() );

                for( size_t side_id = 0; side_id < 2 * level.cmap->dim(); side_id++ )
                {
                    const ElementSide side = elementSideFromId( side_id );
                    const auto adjacent = adjacentElementSide( patch_cmap, local_leaf, side );
                    if( not adjacent.has_value() ) continue;
                    if( not leaves.contains( { patch_id, adjacent->first } ) ) continue;

                    const topology::Cell first_level_cell = level_leaf;
                    const topology::Cell second_level_cell(
                        level.cmap->toGlobalDart( patch_id, adjacent->first.dart() ), level.cmap->dim() );
                    const auto canonical_pair =
                        std::minmax( first_level_cell, second_level_cell );
                    if( not seen_pairs.insert( canonical_pair ).second ) continue;

                    TraceMeshInterface iface{
                        TraceMeshEntityType::Interior,
                        TraceMeshSide{ level.to_space_cell( first_level_cell ), side, patch_id, level.level },
                        TraceMeshSide{ level.to_space_cell( second_level_cell ), adjacent->second, patch_id, level.level },
                        std::nullopt };

                    if( mode == InteriorTraceMode::Broken and
                        not traceHasBrokenRows( ss, iface.first, iface.second.value(), row_tol ) )
                        continue;

                    out.push_back( iface );
                }
            }
        }
        return out;
    }

    std::vector<TraceMeshInterface> traceMeshInterfaces( const SplineSpace& ss,
                                                         const bool include_boundary,
                                                         const bool include_patch_interfaces,
                                                         const bool include_interior,
                                                         const InteriorTraceMode interior_mode,
                                                         const double row_tol )
    {
        std::vector<TraceMeshInterface> out;
        const auto append = [&out]( const std::vector<TraceMeshInterface>& more ) {
            out.insert( out.end(), more.begin(), more.end() );
        };

        if( include_boundary ) append( boundaryTraceMeshInterfaces( ss ) );
        if( include_patch_interfaces ) append( patchTraceMeshInterfaces( ss ) );
        if( include_interior ) append( interiorTraceMeshInterfaces( ss, interior_mode, row_tol ) );
        return out;
    }

    TraceInterfaceElement traceInterfaceElement( const SplineSpace& ss,
                                                 const TraceMeshInterface& iface,
                                                 const double row_tol )
    {
        TraceInterfaceElement out{ traceSideData( ss, iface.first.element, iface.first.side, row_tol ), std::nullopt };
        if( iface.second.has_value() )
        {
            out.second = traceSideData( ss, iface.second->element, iface.second->side, row_tol );
        }
        return out;
    }
}
