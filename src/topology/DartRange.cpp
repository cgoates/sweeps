#include <DartRange.hpp>
#include <TPCombinatorialMap.hpp>
#include <HierarchicalTPCombinatorialMap.hpp>

using namespace topology;

std::pair<size_t, Dart> DartRanges::toLocalDart( const Dart& global_dart ) const
{
    const auto it = std::lower_bound(
        mRanges.begin(), mRanges.end(), global_dart.id(), []( const DartRange& range, const Dart::IndexType d ) {
            return range.max() < d;
        } );
    if( it == mRanges.end() ) throw std::runtime_error( "Dart out of range of multi patch cmap" );

    const size_t patch_id = std::distance( mRanges.begin(), it );
    return { patch_id, mRanges.at( patch_id ).toLocalDart( global_dart ) };
}

Dart DartRanges::toGlobalDart( const size_t i, const Dart& local_d ) const
{
    if( i > mRanges.size() ) throw std::runtime_error( "Requested dart range index is greater than highest range index" );
    return mRanges.at( i ).toGlobalDart( local_d );
}

Dart::IndexType DartRanges::maxDartId() const
{
    return mRanges.back().max();
}

namespace topology
{
    template <typename CMAP>
    DartRanges initializeRanges( const std::vector<std::shared_ptr<const CMAP>>& constituents )
    {
        std::vector<DartRange> out;
        out.reserve( constituents.size() );
        Dart::IndexType next_min = 0;
        for( const auto& patch : constituents )
        {
            out.emplace_back( next_min, *patch );
            next_min = out.back().max() + 1;
        }
        return out;
    }
    template DartRanges initializeRanges<TPCombinatorialMap>( const std::vector<std::shared_ptr<const TPCombinatorialMap>>& );
    template DartRanges initializeRanges<HierarchicalTPCombinatorialMap>( const std::vector<std::shared_ptr<const HierarchicalTPCombinatorialMap>>& );
}