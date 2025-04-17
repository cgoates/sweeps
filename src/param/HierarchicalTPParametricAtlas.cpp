#include <HierarchicalTPParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CustomEigen.hpp>

using namespace param;

HierarchicalTPParametricAtlas::HierarchicalTPParametricAtlas(
    const std::shared_ptr<const topology::HierarchicalTPCombinatorialMap>& cmap,
    const std::vector<std::shared_ptr<const TPParametricAtlas>>& refinement_levels )
    : mMap( cmap ), mRefinementLevels( refinement_levels )
{}

const ParentDomain HierarchicalTPParametricAtlas::parentDomain( const topology::Cell& c ) const
{
    if( c.dim() != mMap->dim() )
        throw std::invalid_argument( "HierarchicalTPParametricAtlas::parentDomain only works for elements." );

    const auto [ level, level_d ] = mMap->unrefinedAncestorDart( c.dart() );
    return mRefinementLevels.at( level )->parentDomain( topology::Cell( level_d, c.dim() ) );
}

ParentPoint HierarchicalTPParametricAtlas::parentPoint( const topology::Vertex& v ) const
{
    if( mMap->dim() == 3 )
    {
        // We want to know how many darts away from the boundary we are in each direction.
        const auto [ dart_level, local_d ] = mMap->dartRanges().toLocalDart( v.dart() );
        const auto [elem_level, unrefined_d] = mMap->unrefinedAncestorDart( v.dart() );

        if( dart_level == elem_level ) return mRefinementLevels.at( elem_level )->parentPoint( topology::Vertex( local_d ) );

        const topology::TPCombinatorialMap& dart_level_cmap = *mMap->refinementLevels().at( dart_level );

        const SmallVector<std::shared_ptr<const ParametricAtlas1d>, 3> pa_1ds =
            tensorProductComponentAtlases( *mRefinementLevels.at( dart_level ) );

        // If we fully unflatten using unflatten cell on the vertex, we will have three vertices.
        const std::array<topology::CellOrEndVertex, 3> fully_unflat_verts = [&]() -> std::array<topology::CellOrEndVertex, 3> {
            const auto unflat1 = unflattenCell( dart_level_cmap, topology::Vertex( local_d ) );
            const std::shared_ptr<const topology::TPCombinatorialMap> source_primal =
                std::dynamic_pointer_cast<const topology::TPCombinatorialMap>( dart_level_cmap.sourceCMapPtr() );
            if( source_primal.get() == nullptr ) throw std::runtime_error( "HierarchicalTPParametricAtlas only supports tensor products of one dimensional parametric atlases" );
            const auto unflat2 = unflattenCell( *source_primal, unflat1.first.value() );
            
            return { unflat2.first, unflat2.second, unflat1.second };
        }();

        const size_t ratio = [&](){
            size_t out = 1;
            for( size_t level_ii = elem_level; level_ii < dart_level; level_ii++ )
                out *= mMap->refinementRatios().at( level_ii );
            return out;
        }();
        const Vector6dMax unrefined_lengths = mRefinementLevels.at( elem_level )->parametricLengths( topology::Cell( unrefined_d, mMap->dim() ) );

        Eigen::Vector3d pt = Eigen::Vector3d::Zero();
        BaryCoordIsZeroVec is_zero;
        for( size_t i = 0; i < 3; i++ )
        {
            if( not fully_unflat_verts.at( i ) )
            {
                pt( i ) = 1.0;
                is_zero.push_back( true );
                is_zero.push_back( false );
            }
            else if( fully_unflat_verts.at( i ).value().dart().id() % ratio == 0 )
            {
                pt( i ) = 0.0;
                is_zero.push_back( false );
                is_zero.push_back( true );
            }
            else
            {
                topology::Dart d = fully_unflat_verts.at( i ).value().dart();
                while( d.id() % ratio != 0 )
                {
                    const size_t len_idx = parametricLengthIndexAlongEdge( *pa_1ds.at( i ), topology::Edge( d ) );
                    pt( i ) += pa_1ds.at( i )->parametricLengths( topology::Cell( d, 1 ) )( len_idx );
                    d = phi( pa_1ds.at( i )->cmap(), -1, d ).value();
                }
                pt( i ) /= unrefined_lengths( i );

                is_zero.push_back( false );
                is_zero.push_back( false );
            }
        }
        return ParentPoint( cubeDomain( 3 ), pt, is_zero );
    }
    else
    {
        // 1. Get the unrefined dart, get the length along that edge.
        const auto [elem_level, unrefined_d] = mMap->unrefinedAncestorDart( v.dart() );

        const TPParametricAtlas& elem_level_atlas = *mRefinementLevels.at( elem_level );
        const double unrefined_len = [&]() {
            const size_t len_idx = parametricLengthIndexAlongEdge( elem_level_atlas, topology::Edge( unrefined_d ) );
            return elem_level_atlas.parametricLengths( topology::Cell( unrefined_d, mMap->dim() ) )( len_idx );
        }();

        // 2. iterate phi(-1) from v.dart() while the answer edge is along the same parametric boundary. Add up the edge lengths of the darts in this iteration.
        double accum_pre_length = 0;
        bool is_zero = true;
        const auto [v_level, level_local_v_dart] = mMap->dartRanges().toLocalDart( v.dart() );
        const BaryCoordIsZeroVec bdry = parentDomainBoundary( *mRefinementLevels.at( v_level ), topology::Edge( level_local_v_dart ) );

        topology::Dart d = phi( *mMap, -1, v.dart() ).value();
        while ( d != v.dart() )
        {
            const auto [d_level, level_d] = mMap->dartRanges().toLocalDart( d );
            if( parentDomainBoundary( *mRefinementLevels.at( d_level ), topology::Edge( level_d ) ) == bdry )
            {
                const size_t len_idx = parametricLengthIndexAlongEdge( *mRefinementLevels.at( d_level ), topology::Edge( level_d ) );
                accum_pre_length += mRefinementLevels.at( d_level )->parametricLengths( topology::Cell( level_d, mMap->dim() ) )( len_idx );
                is_zero = false;
            }
            else
                break;

            d = phi( *mMap, -1, d ).value();
        }

        // 3. Divide the answer from (2) by the answer from (1).
        const ParentPoint bdry_point = [&](){
            if( bdry.at( 3 ) or bdry.at( 0 ) )
                return ParentPoint( cubeDomain( 1 ), Vector1d( accum_pre_length / unrefined_len ), { false, is_zero } );
            else
                return ParentPoint( cubeDomain( 1 ), Vector1d( 1.0 - accum_pre_length / unrefined_len ), { is_zero, false } );
        }();

        // 4. Put the answer from (3) into the appropriate entry in the parent point,
        //    the rest of which comes from parentPoint( Vertex( unrefined_d ) ).
        //    So that's basically just the opposite of a dropCoordinate operation.

        return liftFromBoundary(
            bdry_point,
            mRefinementLevels.at( v_level )->parentDomain( topology::Cell( level_local_v_dart, mMap->dim() ) ),
            bdry );
    }
}

Vector6dMax HierarchicalTPParametricAtlas::parametricLengths( const topology::Cell& c ) const
{
    if( c.dim() != mMap->dim() )
        throw std::invalid_argument( "HierarchicalTPParametricAtlas::parametricLengths only works for elements." );
    const auto [ level, level_d ] = mMap->unrefinedAncestorDart( c.dart() );
    return mRefinementLevels.at( level )->parametricLengths( topology::Cell( level_d, c.dim() ) );
}