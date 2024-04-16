#pragma once
#include <CombinatorialMap.hpp>
#include <GlobalCellMarker.hpp>
#include <CombinatorialMapMethods.hpp>
#include <map>

namespace topology
{
    class LevelSetCMap : public CombinatorialMap
    {
        public:
        LevelSetCMap( const CombinatorialMap& base,
                      const std::function<double( const Vertex& )>& func,
                      const double value );

        LevelSetCMap( const CombinatorialMap& base,
                      const std::function<double( const Vertex& )>& func,
                      const double value,
                      const Edge& one_intersected_edge );

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;
        virtual Dart::IndexType maxDartId() const override { return mUnderlyingMap.maxDartId(); }
        virtual uint dim() const override { return mUnderlyingMap.dim() - 1; }
        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;
        virtual bool iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const override;
        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        private:
        const CombinatorialMap& mUnderlyingMap;
        GlobalCellMarker mIntersectedEdges;
        std::map<Vertex, double> mIntersectionPositions;
    };
}