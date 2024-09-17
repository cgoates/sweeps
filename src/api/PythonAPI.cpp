#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <KnotVector.hpp>
#include <NavierStokesDiscretization.hpp>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE( splines, m )
{
    m.doc() = "Plugin providing a spline discretization for Navier Stokes problems";
    py::class_<basis::KnotVector>( m, "KnotVector" )
        .def( py::init<const std::vector<double>&, const double>(),
              "Constructs a knot vector with the given values. Takes a parametric tolerance to determine which knots "
              "are repeated knots.",
              "knot_values"_a,
              "param_tol"_a )
        .def( "size", &basis::KnotVector::size, "The number of knot values in the knot vector." )
        .def( "knot", &basis::KnotVector::knot, "Returns the knot at position ii.", "ii"_a )
        .doc() = "A simple knot vector object for defining a B-spline.";

    py::enum_<api::PatchSide>( m, "PatchSide" )
        .value( "S0", api::PatchSide::S0 )
        .value( "S1", api::PatchSide::S1 )
        .value( "T0", api::PatchSide::T0 )
        .value( "T1", api::PatchSide::T1 );

    py::class_<topology::Face>( m, "Element" )
        .def( py::init( []( const topology::Dart::IndexType& id ) { return topology::Face( topology::Dart( id ) ); } ),
              "Constructs an element with the given dart id.",
              "dart"_a )
        .def_property_readonly(
            "dart",
            []( const topology::Face& f ) { return f.dart().id(); },
            "Returns the id of a dart (or half edge) on the element." );

    py::class_<topology::Edge>( m, "Edge" )
        .def( py::init( []( const topology::Dart::IndexType& id ) { return topology::Edge( topology::Dart( id ) ); } ),
              "Constructs an edge with the given dart id.",
              "dart"_a )
        .def_property_readonly(
            "dart",
            []( const topology::Edge& e ) { return e.dart().id(); },
            "Returns the id of a dart (or half edge) on the edge." );

    py::class_<eval::SplineSpaceEvaluator>( m, "SplineSpace" )
        .def(
            "numTotalFunctions",
            []( const eval::SplineSpaceEvaluator& eval ) { return eval.splineSpace().numFunctions(); },
            "The total number of functions in the spline space." )
        .def(
            "numVectorComponents",
            []( const eval::SplineSpaceEvaluator& eval ) { return eval.splineSpace().numVectorComponents(); },
            "The number of vector components of the spline space basis." )
        .def(
            "connectivity",
            []( const eval::SplineSpaceEvaluator& eval, const topology::Face& f ) {
                std::vector<basis::FunctionId> conn = eval.splineSpace().connectivity( f );
                std::vector<size_t> out;
                out.reserve( conn.size() );
                std::transform( conn.begin(),
                                conn.end(),
                                std::back_inserter( out ),
                                []( const basis::FunctionId& fid ) { return fid.id(); } );
                return out;
            },
            "Returns the global ids of the basis functions that are nonzero on the given element.",
            "elem"_a )
        .def(
            "numFunctions",
            []( const eval::SplineSpaceEvaluator& eval, const topology::Face& f ) {
                return eval.splineSpace().connectivity( f ).size();
            },
            "Returns the number of nonzero basis functions on the given element.",
            "elem"_a )
        .def( "evaluateBasis",
              &eval::SplineSpaceEvaluator::evaluateBasis,
              "Evaluates the spline space basis at the element and point specified in the most recent calls to "
              "localizeElement and localizePoint on the discretization.  The basis is returned as a matrix with "
              "numFunctions( elem ) rows and numVectorComponents columns, with the nth row corresponding to the nth "
              "global id in the connectivity." )
        .def( "evaluateFirstDerivatives",
              &eval::SplineSpaceEvaluator::evaluateFirstDerivatives,
              "Evaluates the first derivatives of the spine space basis at the element and point specified in the most "
              "recent calls to localizeElement and localizePoint on the discretization. The first derivatives are "
              "returned as a matrix with numFunctions( elem ) rows and numVectorComponents()*2 columns, with the nth "
              "row corresponding to the nth global id in the connectivity, and the columns ordered as the derivative "
              "of all of the vector components with respect to s (one column for each component), followed by the "
              "derivative of all of the vector components with respect to t (another column for each component)." );

    py::class_<api::NavierStokesDiscretization>( m, "NavierStokesDiscretization" )
        .def( py::init<const basis::KnotVector&,
                       const basis::KnotVector&,
                       const size_t,
                       const size_t,
                       const Eigen::Matrix2Xd&>(),
              "Create a B-spline discretization for Navier Stokes problems, consisting of an H1 spline space with the "
              "given knot vectors, an HDiv spline space taking the H1 space as its primal basis, and an L2 spline "
              "space, which has reduced degree from H1 in both directions. The geometry is defined by the H1 space and "
              "the provided control points."
              "knot_vec_s"_a,
              "knot_vec_t"_a,
              "degree_s"_a,
              "degree_t"_a,
              "control_points"_a )
        .def_readonly( "H1", &api::NavierStokesDiscretization::H1, "The H1 spline space." )
        .def_readonly(
            "HDIV", &api::NavierStokesDiscretization::HDIV, "The HDiv spline space, with two vector components." )
        .def_readonly( "L2", &api::NavierStokesDiscretization::L2, "The L2 spline space." )
        .def_readonly( "control_points",
                       &api::NavierStokesDiscretization::cpts,
                       "The control points which define the geometry along with the H1 space." )
        .def(
            "elements",
            []( const api::NavierStokesDiscretization& nsd ) {
                std::vector<topology::Face> out;
                out.reserve( cellCount( nsd.H1_ss.basisComplex().parametricAtlas().cmap(), 2 ) );
                iterateCellsWhile(
                    nsd.H1_ss.basisComplex().parametricAtlas().cmap(), 2, [&out]( const topology::Face& f ) {
                        out.push_back( f );
                        return true;
                    } );
                return out;
            },
            "A list of all the elements in the discretization." )
        .def(
            "boundaryEdges",
            []( const api::NavierStokesDiscretization& nsd ) {
                std::vector<topology::Edge> out;
                out.reserve( cellCount( nsd.cmap_bdry, 1 ) );
                iterateCellsWhile( nsd.cmap_bdry, 1, [&out]( const topology::Edge& e ) {
                    out.push_back( e );
                    return true;
                } );
                return out;
            },
            "A list of all edges on the boundary of the discretization" )
        .def(
            "boundaryEdges",
            []( const api::NavierStokesDiscretization& nsd, const api::PatchSide& side ) {
                const topology::TPCombinatorialMap& tp_map = static_cast<const topology::TPCombinatorialMap&>(
                    nsd.H1.splineSpace().basisComplex().parametricAtlas().cmap() );
                std::vector<topology::Edge> out;
                switch( side )
                {
                    case api::PatchSide::S0:
                    {
                        out.reserve( cellCount( tp_map.lineCMap(), 1 ) );
                        const topology::Dart source_dart( 0 );
                        iterateDartsWhile( tp_map.lineCMap(), [&]( const topology::Dart& line_dart ) {
                            out.push_back(
                                tp_map.flatten( source_dart, line_dart, topology::TPCombinatorialMap::TPDartPos::DartPos3 ) );
                            return true;
                        } );
                        break;
                    }
                    case api::PatchSide::S1:
                    {
                        out.reserve( cellCount( tp_map.lineCMap(), 1 ) );
                        const topology::Dart source_dart( tp_map.sourceCMap().maxDartId() );
                        iterateDartsWhile( tp_map.lineCMap(), [&]( const topology::Dart& line_dart ) {
                            out.push_back(
                                tp_map.flatten( source_dart, line_dart, topology::TPCombinatorialMap::TPDartPos::DartPos1 ) );
                            return true;
                        } );
                        break;
                    }
                    case api::PatchSide::T0:
                    {
                        out.reserve( cellCount( tp_map.sourceCMap(), 1 ) );
                        const topology::Dart line_dart( 0 );
                        iterateDartsWhile( tp_map.sourceCMap(), [&]( const topology::Dart& source_dart ) {
                            out.push_back(
                                tp_map.flatten( source_dart, line_dart, topology::TPCombinatorialMap::TPDartPos::DartPos0 ) );
                            return true;
                        } );
                        break;
                    }
                    case api::PatchSide::T1:
                    {
                        out.reserve( cellCount( tp_map.sourceCMap(), 1 ) );
                        const topology::Dart line_dart( tp_map.sourceCMap().maxDartId() );
                        iterateDartsWhile( tp_map.sourceCMap(), [&]( const topology::Dart& source_dart ) {
                            out.push_back(
                                tp_map.flatten( source_dart, line_dart, topology::TPCombinatorialMap::TPDartPos::DartPos2 ) );
                            return true;
                        } );
                        break;
                    }
                }
                return out;
            },
            "A list of all edges on a given side of the discretization spline patch",
            "side"_a )
        .def(
            "boundaryPerpendicularHDivFuncs",
            []( const api::NavierStokesDiscretization& nsd, const api::PatchSide& side ) {
                const basis::DivConfTPSplineSpace& tp_ss = nsd.HDIV_ss;
                const auto& component_bases = tp_ss.scalarTPBases();
                const bool is_S = side == api::PatchSide::S0 or side == api::PatchSide::S1;

                const util::IndexVec lengths = [&]() {
                    const basis::TPSplineSpace& comp = is_S ? *component_bases.at( 0 ) : *component_bases.at( 1 );
                    const auto component_comps = tensorProductComponentSplines( comp );
                    util::IndexVec out;
                    for( const auto& comp : component_comps ) out.push_back( comp->numFunctions() );
                    return out;
                }();

                const auto iter_dir =
                    [&lengths]( const api::PatchSide& side ) -> SmallVector<std::variant<bool, size_t>, 3> {
                    switch( side )
                    {
                        case api::PatchSide::S0: return { size_t{ 0 }, true };
                        case api::PatchSide::S1: return { lengths.at( 0 ) - 1, true };
                        case api::PatchSide::T0: return { true, size_t{ 0 } };
                        case api::PatchSide::T1: return { true, lengths.at( 1 ) - 1 };
                    }
                }( side );

                std::vector<size_t> out;
                out.reserve( lengths.at( is_S ? 0 : 1 ) );

                const size_t offset = is_S ? 0 : component_bases.at( 0 )->numFunctions();

                util::iterateTensorProduct( lengths, { 0, 1 }, iter_dir, [&]( const util::IndexVec& iv ) {
                    out.emplace_back( util::flatten( iv, lengths ) + offset );
                } );

                return out;
            },
            "A list of all HDIV functions which are nonzero on and perpendicular to a given side of the discretization "
            "spline patch",
            "side"_a )
        .def(
            "localizeElement",
            []( api::NavierStokesDiscretization& nsd, const topology::Face& elem ) {
                nsd.H1.localizeElement( elem );
                nsd.HDIV.localizeElement( elem );
                nsd.L2.localizeElement( elem );
            },
            "Sets all future evaluations of the spline spaces to be on element elem, until localizeElement is called "
            "again.",
            "elem"_a )
        .def(
            "localizePoint",
            []( api::NavierStokesDiscretization& nsd, const Eigen::Vector2d& pt ) {
                const param::ParentPoint ppt( param::cubeDomain( 2 ), pt, { false, false, false, false } );
                nsd.H1.localizePoint( ppt );
                nsd.HDIV.localizePoint( ppt );
                nsd.L2.localizePoint( ppt );
            },
            "Sets all future evaluations of the spline spaces to be on point pt, until localizePoint is called again.",
            "pt"_a )
        .def(
            "mapping",
            []( const api::NavierStokesDiscretization& nsd ) {
                return nsd.H1.evaluateManifold( nsd.cpts );
            },
            "Evaluates the spatial position of the spline geometry at the parametric position from the latest calls to "
            "localizeElement and localizePoint." )
        .def(
            "jacobianDeterminant",
            []( const api::NavierStokesDiscretization& nsd ) {
                return nsd.H1.evaluateJacobian( nsd.cpts ).determinant();
            },
            "Evaluates the Jacobian determinant of the spline geometry at the parametric position from the latest "
            "calls to localizeElement and LocalizePoint." )
        .def(
            "piolaTransformedHDIVBasis",
            []( const api::NavierStokesDiscretization& nsd ) {
                return piolaTransformedVectorBasis( nsd.HDIV, nsd.H1, nsd.cpts );
            },
            "Evaluates the Piola transformed HDiv basis at the parametric position from the latest calls to "
            "localizeElement and localizePoint. The basis is returned as a matrix with HDIV.numFunctions( elem ) rows "
            "and two columns, with the ith row corresponding to the ith global id in the connectivity, and the jth "
            "column corresponding to the jth component of the vector valued basis." )
        .def(
            "piolaTransformedHDIVFirstDerivatives",
            []( const api::NavierStokesDiscretization& nsd ) {
                return piolaTransformedVectorFirstDerivatives( nsd.HDIV, nsd.H1, nsd.cpts );
            },
            "Evaluates the Piola transformed first derivatives of the HDiv basis at the parametric position from the "
            "latest calls to localizeElement and localizePoint. The basis is returned as a matrix with "
            "HDIV.numFunctions( elem ) rows and four columns, with the ith row corresponding to the ith global id in "
            "the connectivity, and the columns ordered as dv_x/ds, dv_y/ds, dv_x/dt, dv_y/dt." )
        .def(
            "piolaTransformedL2",
            []( const api::NavierStokesDiscretization& nsd ) {
                return piolaTransformedBivectorBasis( nsd.L2, nsd.H1, nsd.cpts );
            },
            "Evaluates the Piola transformed L2 basis at the parametric position from the latest calls to "
            "localizeElement and localizePoint. The basis is returned as a matrix with L2.numFunctions( elem ) rows "
            "and one column, with the ith row corresponding to the ith global id in the connectivity." )
        .def(
            "piolaTransformedL2FirstDerivatives",
            []( const api::NavierStokesDiscretization& nsd ) {
                return piolaTransformedBivectorFirstDerivatives( nsd.L2, nsd.H1, nsd.cpts );
            },
            "Evaluates the piola transformed derivatives of the L2 basis at the parametric position from the latest "
            "calls to localizeElement and localizePoint. The basis is returned as a matrix with L2.numFunctions( elem "
            ") rows and two columns, with the ith row corresponding to the ith global id in the connectivity, and the "
            "columns ordered as dp/ds, dp/dt." );

    m.def(
        "grevillePoints",
        []( const basis::KnotVector& kv_s, const basis::KnotVector& kv_t, const size_t degree_s, const size_t degree_t )
            -> Eigen::Matrix2Xd {
            return util::tensorProduct( { grevillePoints( kv_s, degree_s ), grevillePoints( kv_t, degree_t ) } ).transpose();
        },
        "Returns the greville points of the 2d B-spline patch defined by the given knot vectors and degrees.  These "
        "points, when used as control points, create a linear spatial spline geometry corresponding exactly to the "
        "parametic domain."
        "kv_s"_a,
        "kv_t"_a,
        "degree_s"_a,
        "degree_t"_a );
}