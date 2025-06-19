#include <SweepAPIMethods.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <SweepInput.hpp>
#include <MeshInput.hpp>

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

    py::class_<api::HexMesh>( m, "HexMesh" )
        .def_readwrite( "points", &api::HexMesh::points )
        .def_readwrite( "hexes", &api::HexMesh::hexes )
        .doc() = "A simple hex mesh, with a list of points, and a list of hexes.  The vertices of the hexes are "
                 "ordered by their coordinates in the hex-local coordinates system as (0,0,0), (1,0,0), (0,1,0), "
                 "(1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1).";

    py::class_<api::QuadMesh>( m, "QuadMesh" )
        .def_readwrite( "points", &api::QuadMesh::points )
        .def_readwrite( "quads", &api::QuadMesh::quads )
        .doc() = "A simple quad mesh, with a list of points, and a list of quads.  The vertices of the quads are "
                 "ordered by their coordinates in the local coordinates system as (0,0), (1,0), (1,1), (0,1).";

    py::class_<api::TriMesh>( m, "TriMesh" )
        .def_readwrite( "points", &api::TriMesh::points )
        .def_readwrite( "tris", &api::TriMesh::tris )
        .doc() = "A simple triangle mesh, with a list of points, and a list of triangles.";

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

    m.def(
        "fitSinglePatchHexMeshToSweep",
        &api::fitSinglePatchHexMeshToSweep,
        "Fits a single patch hex mesh to a sweep.",
        "sweep"_a,
        "n_elems_st"_a,
        "u_values"_a,
        py::kw_only(),
        "debug"_a = false
    );

    m.def(
        "fitFivePatchHexMeshToSweep",
        &api::fitFivePatchHexMeshToSweep,
        "Fits a five patch hex mesh to a sweep.",
        "sweep"_a,
        "n_elems_st"_a,
        "u_values"_a,
        py::kw_only(),
        "debug"_a = false
    );

    m.def(
        "fitHexMeshToSweep",
        &api::fitHexMeshToSweep,
        "Fits a hex mesh to a sweep, starting from a user-specified quad mesh on the source surface.",
        "sweep"_a,
        "quad_mesh"_a,
        "u_values"_a,
        py::kw_only(),
        "debug"_a = false
    );

    m.def(
        "loadQuadMeshFromObjFile",
        []( const std::string& filename ){
            const auto [quads, points] = io::loadOBJFile( filename );
            std::vector<std::array<VertexId::Type, 4>> quads_out;
            quads_out.reserve( quads.size() );
            for( const auto& q : quads )
            {
                quads_out.push_back( { q[0].id(), q[1].id(), q[2].id(), q[3].id() } );
            }
            return api::QuadMesh{ points, quads_out };
        },
        "filename"_a,
        "Loads a quad mesh from an OBJ file. The OBJ file should contain quads, not triangles."
    );

    m.def(
        "baseOfSweep",
        &api::baseOfSweep,
        "Computes the base surface triangle mesh of a sweep.",
        "sweep"_a
    );
}