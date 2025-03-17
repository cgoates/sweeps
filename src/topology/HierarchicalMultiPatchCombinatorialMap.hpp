#pragma once
#include <MultiPatchCombinatorialMap.hpp>
#include <HierarchicalTPCombinatorialMap.hpp>
#include <DartRange.hpp>
#include <map>

namespace topology
{
    /// @brief A hierarchical combinatorial map over multi-patch topologies.
    ///
    class HierarchicalMultiPatchCombinatorialMap : public CombinatorialMap
    {
        public:
        HierarchicalMultiPatchCombinatorialMap( const std::vector<std::shared_ptr<const MultiPatchCombinatorialMap>>& refinement_levels,
                                                const std::vector<std::vector<Cell>>& leaf_elements );
        virtual ~HierarchicalMultiPatchCombinatorialMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override;

        virtual uint dim() const override;

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;

        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint ) const override;

        std::pair<size_t, Dart> unrefinedAncestorDart( const Dart& leaf_d ) const;

        const std::vector<std::shared_ptr<const MultiPatchCombinatorialMap>>& refinementLevels() const { return mRefinementLevels; }

        const std::vector<std::shared_ptr<const HierarchicalTPCombinatorialMap>>& constituents() const { return mConstituents; }

        size_t numLevels() const { return mRefinementLevels.size(); }

        const DartRanges& dartRanges() const { return mRanges; }

        bool iterateChildren( const Cell& local_cell,
                              const size_t cell_level,
                              const std::function<bool( const Cell& )>& callback ) const;

        private:
        const std::vector<std::shared_ptr<const MultiPatchCombinatorialMap>> mRefinementLevels;
        std::vector<std::shared_ptr<const HierarchicalTPCombinatorialMap>> mConstituents;
        const DartRanges mRanges;
    };

    std::vector<std::vector<Cell>> leafElements( const HierarchicalMultiPatchCombinatorialMap& cmap );
}; // namespace topology