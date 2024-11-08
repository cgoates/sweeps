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

    void outputLevelSetsAndTraces( const Sweep& sweep,
                                   const std::vector<double>& level_set_values,
                                   const std::vector<Eigen::Vector2d>& trace_points,
                                   const std::string& output_prefix );
}