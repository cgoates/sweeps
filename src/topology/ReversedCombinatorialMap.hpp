#pragma once
#include <CombinatorialMap.hpp>
#include <set>
#include <VertexPositionsFunc.hpp>

namespace topology
{
    /// @brief A combinatorial map representing the boundary of another combinatorial map.
    ///
    class ReversedCombinatorialMap : public CombinatorialMap
    {
        public:
        explicit ReversedCombinatorialMap( const CombinatorialMap& map );
        virtual ~ReversedCombinatorialMap() = default;
        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;
        virtual Dart::IndexType maxDartId() const override;
        virtual uint dim() const override;
        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;
        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;
        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        topology::Cell toUnderlyingCell( const topology::Cell& c ) const;
        topology::Cell fromUnderlyingCell( const topology::Cell& c ) const;

        private:
        const CombinatorialMap& mUnderlyingMap;
    };

    VertexPositionsFunc reversedVertexPositions( const ReversedCombinatorialMap& bdry,
                                                 const VertexPositionsFunc& underlying_positions );
}; // namespace topology