#pragma once
#include <SweepInput.hpp>
#include <optional>

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

template <unsigned int DIM> struct Triangle;
template <unsigned int DIM> struct Segment;
class Normal;

template<unsigned int DIM>
struct Ray
{
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> start_pos;
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> dir;
};

std::optional<Eigen::Vector3d> intersectionOf( const Ray<3>& ray,
                                               const Triangle<3>& tri,
                                               std::optional<const Eigen::Vector3d> maybe_normal = {} );

using TracePoint = std::pair<topology::Face, Eigen::Vector3d>;
std::optional<TracePoint> traceRayOnTet( const topology::TetMeshCombinatorialMap& map,
                                         const topology::Volume& v,
                                         const Ray<3>& ray,
                                         const std::vector<Normal>& normals );

SimplicialComplex traceField( const topology::TetMeshCombinatorialMap& map,
                              const topology::Cell& start_cell,
                              const Eigen::Vector3d& start_point,
                              const Eigen::Matrix3Xd& field,
                              const std::vector<Normal>& normals,
                              const bool debug_output = false );

std::optional<Eigen::Vector2d> intersectionOf( const Ray<2>& ray,
                                               const Segment<2>& line );

std::optional<std::pair<bool, double>> traceGradientOnTri( const Triangle<3>& tri3d,
                                                           const double edge_barycentric_coord,
                                                           const Eigen::Ref<const Eigen::Vector3d> field_values );

using VertexPositionsFunc = std::function<const Eigen::Vector3d&( const topology::Vertex& )>;
std::optional<std::pair<topology::Edge, double>> traceGradientOnTri( const topology::CombinatorialMap& map,
                                                                     const VertexPositionsFunc& positions,
                                                                     const topology::Face& f,
                                                                     const double edge_barycentric_coord,
                                                                     const Eigen::VectorXd& field_values );

SimplicialComplex traceBoundaryField( const topology::CombinatorialMap& map,
                                      const topology::Edge& e,
                                      const double& start_point,
                                      const Eigen::VectorXd& field,
                                      const std::function<const Eigen::Vector3d&( const topology::Vertex& )>& positions,
                                      const bool debug_output,
                                      const std::function<void( const topology::Face& )>& face_callback = []( const auto& ) {} );