#pragma once
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/core/types/cell_marker.h>
#include <SweepInput.hpp>

struct Triangle;
class Normal;

struct Ray
{
    const Eigen::Ref<const Eigen::Vector3d> start_pos;
    const Eigen::Ref<const Eigen::Vector3d> dir;
};
Eigen::MatrixX3d gradients( const cgogn::CMap3& map,
                            const Eigen::VectorXd& field_values,
                            const std::vector<Normal>& normals );

std::optional<Eigen::Vector3d> intersectionOf( const Ray& ray,
                                               const Triangle& tri,
                                               std::optional<const Eigen::Vector3d> maybe_normal = {} );

using TracePoint = std::pair<cgogn::CMap3::Face, Eigen::Vector3d>;
std::optional<TracePoint> traceRayOnTet( const cgogn::CMap3& map,
                                         const cgogn::CMap3::Volume& v,
                                         const Ray& ray,
                                         const std::vector<Normal>& normals );

SimplicialComplex traceField( const cgogn::CMap3& map,
                              const cgogn::CMap3::Face& f,
                              const Eigen::Vector3d& start_point,
                              const Eigen::MatrixX3d& field,
                              const std::vector<Normal>& normals,
                              const bool debug_output = false );