#pragma once
#include <map>
#include <vector>
#include <optional>
#include <SmallVector.hpp>
#include <CombinatorialMap.hpp>

namespace topology
{
    /// @brief A helper class for explicitly constructing a combinatorial map.
    /// Has not been optimized, should only be used for testing situations.
    class CustomCombinatorialMap : public CombinatorialMap
    {
        public:
        CustomCombinatorialMap( const size_t n_darts,
                                const uint dim,
                                const SmallVector<std::vector<std::optional<Dart::IndexType>>, 3>& phis,
                                const std::map<size_t, std::vector<size_t>>& cell_ids );

        virtual ~CustomCombinatorialMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;
        virtual Dart::IndexType maxDartId() const override;
        virtual uint dim() const override;
        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;
        virtual bool iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const override;
        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;
        virtual std::optional<size_t> cellCount( const uint cell_dim ) const override;

        private:
        std::vector<std::optional<Dart::IndexType>> mPhi_1s;
        std::vector<std::optional<Dart::IndexType>> mPhi1s;
        std::vector<std::optional<Dart::IndexType>> mPhi2s;
        std::vector<std::optional<Dart::IndexType>> mPhi3s;

        std::map<size_t, std::vector<size_t>> mCellIds;

        size_t mDim;
    };

} // namespace topology