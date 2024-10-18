#pragma once
#include <CombinatorialMap.hpp>
#include <map>
#include <VertexPositionsFunc.hpp>

namespace topology
{
    class DelaunayTriangulation : public CombinatorialMap
    {
        public:
        DelaunayTriangulation( const std::shared_ptr<const CombinatorialMap>& base, const VertexPositionsFunc& vert_positions );

        virtual ~DelaunayTriangulation() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override { return mMaxDartId; }

        virtual uint dim() const override { return 2; }

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;
        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint cell_dim ) const override;

        const CombinatorialMap& baseMap() const { return *mBaseMap; }

        private:
        const std::shared_ptr<const CombinatorialMap> mBaseMap;
        std::map<Dart, Dart> mAlteredPhi1s;
        std::map<Dart, Dart> mAlteredPhi_1s;
        const Dart::IndexType mLowerBound;
        Dart::IndexType mMaxDartId;
    };

    VertexPositionsFunc delaunayTriangulationVertexPositions(
        const DelaunayTriangulation& dtri, const VertexPositionsFunc& underlying_positions );
} // namespace topology