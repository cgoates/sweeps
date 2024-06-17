#pragma once
#include <CombinatorialMap.hpp>

class VertexId;

namespace topology
{
    /// @brief A one-dimensioncal combinatorial map with a specified number of cells
    /// that may or may not be periodic.
    class CombinatorialMap1d : public CombinatorialMap
    {
        public:
        CombinatorialMap1d( const Dart::IndexType n_cells, const bool periodic = false );
        virtual ~CombinatorialMap1d() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override;

        virtual uint dim() const override;

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;

        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint ) const override;

        private:
        Dart::IndexType mNumCells;
        bool mPeriodic;
    };
}; // namespace topology