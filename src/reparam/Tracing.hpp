#pragma once
#include <SimplicialComplex.hpp>
#include <optional>
#include <VertexPositionsFunc.hpp>

namespace topology
{
    class CombinatorialMap;
    class TetMeshCombinatorialMap;
    class Volume;
    class Face;
    class Edge;
    class Vertex;
    class Cell;
}

template <int DIM> struct Triangle;
template <int DIM> struct Segment;
class Normal;

template <int DIM> struct Ray
{
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> start_pos;
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> dir;
};

namespace reparam
{
    class Trace
    {
        public:
        SimplicialComplex mComplex;
        std::vector<double> mHarmonicFuncValues;
        std::vector<topology::Cell> mBaseCells;
    };

    std::optional<Eigen::Vector3d> intersectionOf( const Ray<3>& ray,
                                                   const Triangle<3>& tri,
                                                   std::optional<const Eigen::Vector3d> maybe_normal = {} );

    using TracePoint = std::pair<topology::Face, Eigen::Vector3d>;
    std::optional<TracePoint> traceRayOnTet( const topology::TetMeshCombinatorialMap& map,
                                             const topology::Cell& start_cell,
                                             const Ray<3>& ray,
                                             const std::vector<Normal>& normals );

    SimplicialComplex traceField( const topology::TetMeshCombinatorialMap& map,
                                  const topology::Cell& start_cell,
                                  const Eigen::Vector3d& start_point,
                                  const Eigen::Matrix3Xd& field,
                                  const std::vector<Normal>& normals,
                                  const bool debug_output = false );

    std::optional<Eigen::Vector2d> intersectionOf( const Ray<2>& ray, const Segment<2>& line );

    std::optional<std::pair<topology::Edge, double>> traceGradientOnTri( const topology::CombinatorialMap& map,
                                                                         const VertexPositionsFunc& positions,
                                                                         const topology::Cell& start_cell,
                                                                         const double edge_barycentric_coord,
                                                                         const Eigen::VectorXd& field_values );

    Trace traceBoundaryField( const topology::CombinatorialMap& map,
                              const topology::Cell& start_cell,
                              const double& start_point,
                              const Eigen::VectorXd& field,
                              const VertexPositionsFunc& positions,
                              const bool debug_output );
} // namespace reparam