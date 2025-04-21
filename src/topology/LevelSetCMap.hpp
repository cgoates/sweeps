#pragma once
#include <CombinatorialMap.hpp>
#include <GlobalCellMarker.hpp>
#include <CombinatorialMapMethods.hpp>
#include <map>
#include <Eigen/Core>

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

        virtual ~LevelSetCMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;
        virtual Dart::IndexType maxDartId() const override { return mUnderlyingMap.maxDartId(); }
        virtual uint dim() const override { return mUnderlyingMap.dim() - 1; }
        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;
        virtual bool iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const override;
        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        double intersectionPosition( const topology::Vertex& v ) const;
        Cell underlyingCell( const topology::Cell& c ) const;
        const CombinatorialMap& underlyingMap() const { return mUnderlyingMap; }

        private:
        const CombinatorialMap& mUnderlyingMap;
        GlobalCellMarker mIntersectedEdges;
        std::map<Vertex, double> mIntersectionPositions;
        std::map<size_t, size_t> mVertexIds; // Remove if we no longer need contiguous zero-based ids
    };

    std::function<Eigen::Vector3d( const Vertex& )> levelSetVertexPositions(
        const LevelSetCMap& level, const std::function<Eigen::Vector3d( const Vertex& )>& underlying_positions );
}