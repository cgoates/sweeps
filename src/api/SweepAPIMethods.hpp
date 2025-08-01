#pragma once
#include <vector>
#include <optional>
#include <string>
#include <set>
#include <SimplicialComplex.hpp>

struct SweepInput;

namespace api
{
    struct Sweep
    {
        SimplicialComplex mesh;
        std::set<VertexId::Type> source;
        std::set<VertexId::Type> target;
    };

    struct HexMesh
    {
        std::vector<Eigen::Vector3d> points;
        std::vector<std::array<VertexId::Type, 8>> hexes;
    };

    struct QuadMesh
    {
        std::vector<Eigen::Vector3d> points;
        std::vector<std::array<VertexId::Type, 4>> quads;
    };

    struct TriMesh
    {
        std::vector<Eigen::Vector3d> points;
        std::vector<std::array<VertexId::Type, 3>> tris;
    };

    void outputLevelSetsAndTraces( const Sweep& sweep,
                                   const std::vector<double>& level_set_values,
                                   const std::vector<Eigen::Vector2d>& trace_points,
                                   const std::string& output_prefix );

    HexMesh fitSinglePatchHexMeshToSweep( const api::Sweep& sweep, const size_t n_elems_st, const std::vector<double>& u_values, const bool debug = false );

    HexMesh fitFivePatchHexMeshToSweep( const api::Sweep& sweep, const size_t n_elems_st, const std::vector<double>& u_values, const bool debug = false );

    HexMesh fitHexMeshToSweep( const api::Sweep& sweep, const api::QuadMesh& quad_mesh, const std::vector<double>& u_values, const bool debug = false );

    TriMesh baseOfSweep( const Sweep& sweep );
}