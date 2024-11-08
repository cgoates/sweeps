#include <SweepAPIMethods.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <SweepInput.hpp>
#include <AbaqusInput.hpp>

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE( sweeps, m )
{
    m.doc() = "Plugin providing functions for performing swept parameterization";
    py::class_<Simplex>( m, "Tet" )
        .def_property_readonly(
            "vertices",
            []( const Simplex& s ){
                return std::vector<VertexId::Type>({ s.vertex( 0 ).id(), s.vertex( 1 ).id(), s.vertex( 2 ).id(), s.vertex( 3 ).id() } );
            } )
        .doc() = "A tet, part of a tet mesh. Has four vertices, specified by their id in the tet mesh.";

    py::class_<SimplicialComplex>( m, "TetMesh" )
        .def_readwrite(
            "points",
            &SimplicialComplex::points )
        .def_readwrite(
            "tets",
            &SimplicialComplex::simplices )
        .doc() = "A simple tet mesh, with a list of points, and a list of tets.";

    py::class_<api::Sweep>( m, "Sweep" )
        .def_readwrite(
            "mesh",
            &api::Sweep::mesh )
        .def_readwrite(
            "source",
            &api::Sweep::source )
        .def_readwrite(
            "target",
            &api::Sweep::target )
        .doc() = "A simple sweep specification object, including a tet mesh, and specifications of the sweep source and target as sets of vertex ids.";

    m.def(
        "writeParameterizationToFile",
        &api::outputLevelSetsAndTraces,
        "Computes the level sets of the sweep parameterization and traces the parameterization, outputting the results to .vtu files.",
        "sweep"_a,
        "level_set_values"_a,
        "trace_points"_a,
        "output_prefix"_a );

    m.def(
        "loadFromFile",
        []( const std::string& filename, const std::string& source_set, const std::string& target_set ){
            const SweepInput si = io::loadINPFile( filename, source_set, target_set );
            api::Sweep out{ si.mesh, {}, {} };
            for( size_t i = 0; i < si.mesh.points.size(); i++ )
            {
                if( si.zero_bcs.at( i ) ) out.source.insert( i );
                if( si.one_bcs.at( i ) ) out.target.insert( i );
            }

            return out;
        },
        "filename"_a,
        "source_set_label"_a,
        "target_set_label"_a );
}