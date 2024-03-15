#pragma once
#include <Cell.hpp>
#include <functional>

class VertexId;

namespace topology
{
    class CombinatorialMap
    {
        public:
        virtual std::optional<Dart> phi( const int i, const Dart& d ) const = 0;
        virtual Dart::IndexType maxDartId() const = 0;
        virtual uint dim() const = 0;
        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const = 0;
        virtual bool iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const = 0;

        virtual VertexId vertexId( const Vertex& v ) const = 0;
    };
}; // namespace topology