#pragma once
#include <CombinatorialMap.hpp>
#include <set>
#include <map>
#include <VertexPositionsFunc.hpp>

namespace topology
{
    /// @brief A combinatorial map layer that takes another cmap and cuts it along a specified set of (d-1)-cells.
    /// All darts remain the same, but phi operations change near the cuts, and indexing of vertices is augmented accordingly.
    class CutCombinatorialMap : public CombinatorialMap
    {
        public:
        CutCombinatorialMap( const CombinatorialMap& base, const std::set<topology::Cell>& interfaces_to_disconnect );

        virtual ~CutCombinatorialMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override { return mBaseMap.maxDartId(); }

        virtual uint dim() const override { return mBaseMap.dim(); }

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;

        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint cell_dim ) const override;

        const CombinatorialMap& baseMap() const { return mBaseMap; }

        private:
        const CombinatorialMap& mBaseMap;
        std::set<Dart> mNoPhiDimDarts;
        size_t mNumCuts;
        std::map<Vertex, size_t> mAdditionalVertexIds;
    };
} // namespace topology