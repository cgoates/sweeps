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
    const auto [ level, level_d ] = mMap->unrefinedAncestorDart( c.dart() );
    return mRefinementLevels.at( level )->parentDomain( topology::Cell( level_d, c.dim() ) );
}

ParentPoint HierarchicalTPParametricAtlas::parentPoint( const topology::Vertex& v ) const
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

Vector6dMax HierarchicalTPParametricAtlas::parametricLengths( const topology::Cell& c ) const
{
    const auto [ level, level_d ] = mMap->unrefinedAncestorDart( c.dart() );
    return mRefinementLevels.at( level )->parametricLengths( topology::Cell( level_d, c.dim() ) );
}