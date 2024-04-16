#pragma once
#include <CombinatorialMap.hpp>
#include <GlobalDartMarker.hpp>

class VertexId;

namespace topology
{
    class CombinatorialMapRestriction : public CombinatorialMap
    {
        public:
        /// @param map The map to restrict
        /// @param restriction_func A restriction function called on each element. If true, that element is included.
        CombinatorialMapRestriction( const CombinatorialMap& map,
                                     const std::function<bool( const Cell& )>& restriction_func );

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;
        virtual Dart::IndexType maxDartId() const override { return mUnrestrictedMap.maxDartId(); }
        virtual uint dim() const override { return mUnrestrictedMap.dim(); }
        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;
        virtual bool iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const override;
        virtual VertexId vertexId( const Vertex& v ) const override;
        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        private:
        const CombinatorialMap& mUnrestrictedMap;
        GlobalDartMarker mIncludedDarts;
    };
}