#pragma once
#include <MultiPatchCombinatorialMap.hpp>

namespace topology
{
    class IndexingCellMarker;

    struct MultiPatchDecomposition
    {
        std::vector<std::shared_ptr<const TPCombinatorialMap>> constituents;
        std::map<std::pair<size_t, Dart>, std::pair<size_t, Dart>> connections;
        Dart unstructured_first_corner;
    };

    IndexingCellMarker multiPatchCorners( const CombinatorialMap& unstructured_cmap );
    MultiPatchDecomposition multiPatchDecomposition( const CombinatorialMap& unstructured_cmap );
}