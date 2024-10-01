#pragma once
#include <string>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>

namespace io
{
//    class OBJOutputObject
//    {
//        public:
//        OBJOutputObject( const topology::CombinatorialMap& cmap ) : mCmap( cmap ) {}
//
//        const CombinatorialMap& cmap() const { return mCmap; }
//
//        private:
//        const CombinatorialMap& mCmap;
//    };
//
    /// @brief write out a simplicial mesh to OBJ file format
    /// @param output the output object
    /// @param filename The output filename, including the extension.
    void outputSimplicialMeshToOBJ( const topology::TetMeshCombinatorialMap& to_output, const std::string& filename );

} // namespace io
