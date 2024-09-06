#pragma once
#include <TPCombinatorialMap.hpp>
#include <DartRange.hpp>
#include <map>

namespace topology
{
    /// @brief A hierarchical combinatorial map over tensor product topologies.
    /// Designed to work only with tensor products of 1d cmaps, i.e. not sweeps of unstructured 2d cmaps.
    class HierarchicalTPCombinatorialMap : public CombinatorialMap
    {
        public:
        HierarchicalTPCombinatorialMap( const std::vector<std::shared_ptr<const TPCombinatorialMap>>& refinement_levels,
                                        const std::vector<std::vector<Cell>>& leaf_elements );
        virtual ~HierarchicalTPCombinatorialMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override;

        virtual uint dim() const override;

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;

        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint ) const override;

        std::pair<size_t, Dart> unrefinedAncestorDart( const Dart& leaf_d ) const;

        const DartRanges& dartRanges() const { return mRanges; }

        private:
        bool iterateDartLineage( const Dart& global_d,
                                 const size_t ancestor_or_descendant_level,
                                 const std::function<bool( const Dart& )>& callback ) const;

        bool iterateAncestors( const Dart& global_d,
                               const std::function<bool( const Dart& )>& callback ) const;

        bool iterateLeafDescendants( const Dart& global_d, const std::function<bool( const Dart& )>& callback ) const;

        const std::vector<std::shared_ptr<const TPCombinatorialMap>> mRefinementLevels;
        const DartRanges mRanges;
        std::vector<bool> mLeafDarts;
        std::vector<bool> mUnrefinedDarts;
        std::map<Dart, Dart> mPhiOnes;
        std::map<Dart, Dart> mPhiMinusOnes;
    };
}; // namespace topology