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

        const std::vector<std::shared_ptr<const TPCombinatorialMap>>& refinementLevels() const { return mRefinementLevels; }

        size_t numLevels() const { return mRefinementLevels.size(); }

        bool iterateChildren( const Cell& local_cell,
                              const size_t cell_level,
                              const std::function<bool( const Cell& )>& callback ) const;

        bool iterateLeafDescendants( const Dart& global_d, const std::function<bool( const Dart& )>& callback ) const;

        const std::vector<size_t>& refinementRatios() const { return mRefinementRatios; }

        protected:
        bool iterateDartLineage( const Dart& global_d,
                                 const size_t ancestor_or_descendant_level,
                                 const std::function<bool( const Dart& )>& callback ) const;

        bool iterateAncestors( const Dart& global_d,
                               const std::function<bool( const Dart& )>& callback ) const;

        const std::vector<std::shared_ptr<const TPCombinatorialMap>> mRefinementLevels;
        const DartRanges mRanges;
        std::vector<bool> mLeafDarts;
        std::vector<bool> mUnrefinedDarts;
        std::map<Dart, Dart> mPhiOnes;
        std::map<Dart, Dart> mPhiMinusOnes;
        std::vector<size_t> mRefinementRatios;
    };

    /// @brief A mutable version of HierarchicalTPCombinatorialMap.
    /// This is to be used only for initialization, then call asImmutable() to get an immutable version.
    class MutableHierarchicalTPCombinatorialMap : public HierarchicalTPCombinatorialMap
    {
        public:
        MutableHierarchicalTPCombinatorialMap(
            const std::vector<std::shared_ptr<const TPCombinatorialMap>>& refinement_levels,
            const std::vector<std::vector<Cell>>& leaf_elements )
            : HierarchicalTPCombinatorialMap( refinement_levels, leaf_elements )
        {}

        void setPhi( const Dart& d, const Dart& phi_one )
        {
            mPhiOnes.insert_or_assign( d, phi_one );
            mPhiMinusOnes.insert_or_assign( phi_one, d );
        }
        void setLeaf( const Dart& d, const bool is_leaf )
        {
            mLeafDarts.at( d.id() ) = is_leaf;
        }

        bool isLeaf( const Dart& d ) const
        {
            return mLeafDarts.at( d.id() );
        }

        size_t refinementRatio( const size_t level ) const { return mRefinementRatios.at( level ); }

        /// Expose a protected method from the base class for initialization purposes.
        bool iterateDartLineage( const Dart& global_d,
                                 const size_t ancestor_or_descendant_level,
                                 const std::function<bool( const Dart& )>& callback ) const
        {
            return HierarchicalTPCombinatorialMap::iterateDartLineage( global_d, ancestor_or_descendant_level, callback );
        }

        /// Expose a protected method from the base class for initialization purposes.
        bool iterateAncestors( const Dart& global_d,
                               const std::function<bool( const Dart& )>& callback ) const
        {
            return HierarchicalTPCombinatorialMap::iterateAncestors( global_d, callback );
        }

        HierarchicalTPCombinatorialMap asImmutable() const
        {
            return HierarchicalTPCombinatorialMap( *this );
        }
    };

    bool checkForNoAncestor( const TPCombinatorialMap& tp_map, const Dart& d, const size_t n_darts_per_ancestor );

    std::vector<std::vector<Cell>> leafElements( const HierarchicalTPCombinatorialMap& cmap );
}; // namespace topology