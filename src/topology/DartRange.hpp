#pragma once
#include <Dart.hpp>
#include <vector>
#include <CombinatorialMap.hpp>

namespace topology
{
    class TPCombinatorialMap;

    class DartRange
    {
        public:
        DartRange( const Dart::IndexType min_id, const Dart::IndexType max_id ) : mMin( min_id ), mMax( max_id ) {}
        DartRange( const Dart::IndexType min_id, const CombinatorialMap& cmap )
            : DartRange( min_id, min_id + cmap.maxDartId() )
        {}

        bool contains( const Dart& d ) const { return mMin <= d.id() and mMax >= d.id(); }

        Dart toGlobalDart( const Dart& local_d ) const { return Dart( mMin + local_d.id() ); }

        Dart toLocalDart( const Dart& global_d ) const { return Dart( global_d.id() - mMin ); }

        Dart::IndexType min() const { return mMin; }
        Dart::IndexType max() const { return mMax; }

        private:
        const Dart::IndexType mMin;
        const Dart::IndexType mMax;
    };

    class DartRanges
    {
        public:
        DartRanges( const std::vector<DartRange>& ranges ) : mRanges( ranges ) {}

        std::pair<size_t, Dart> toLocalDart( const Dart& global_d ) const;

        Dart toGlobalDart( const size_t i, const Dart& local_d ) const;

        Dart::IndexType maxDartId() const;

        private:
        std::vector<DartRange> mRanges;
    };

    DartRanges initializeRanges( const std::vector<std::shared_ptr<const TPCombinatorialMap>>& constituents );
}