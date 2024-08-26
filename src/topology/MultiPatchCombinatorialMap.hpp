#pragma once
#include <CombinatorialMap.hpp>
#include <TPCombinatorialMap.hpp>
#include <map>

namespace topology
{
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

    class MultiPatchCombinatorialMap : public CombinatorialMap
    {
        public:
        MultiPatchCombinatorialMap( const std::vector<std::shared_ptr<const TPCombinatorialMap>>& constituents,
                                    const std::map<std::pair<size_t, Dart>, std::pair<size_t, Dart>>& connections );

        virtual ~MultiPatchCombinatorialMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override;

        virtual uint dim() const override;

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;

        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint ) const override;

        struct ConstituentSide
        {
            size_t constituent_id;
            size_t side_id;

            bool operator<( const ConstituentSide& o ) const;
        };

        enum class TPPermutation : Dart::IndexType
        {
            ZeroToZero,
            ZeroToOne,
            ZeroToTwo,
            ZeroToThree,
            Flip1d
        };

        std::pair<size_t, Dart> toLocalDart( const Dart& global_dart ) const;
        Dart toGlobalDart( const size_t patch_id, const Dart& local_dart ) const;

        private:
        std::vector<std::shared_ptr<const TPCombinatorialMap>> mSubMaps;
        std::map<ConstituentSide, std::pair<TPPermutation, ConstituentSide>> mInterMapConnections;
        std::vector<DartRange> mRanges;
    };
} // namespace topology