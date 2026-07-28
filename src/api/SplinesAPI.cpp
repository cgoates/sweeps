#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <KnotVector.hpp>
#include <NavierStokesDiscretization.hpp>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>
#include <Eigen/LU>

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
        .def( "__repr__", []( const basis::KnotVector& kv ) {
            std::stringstream ss;
            ss << kv;
            return ss.str();
        } )
        .doc() = "A simple knot vector object for defining a B-spline.";

    py::enum_<api::PatchSide>( m, "PatchSide" )
        .value( "S0", api::PatchSide::S0 )
        .value( "S1", api::PatchSide::S1 )
        .value( "T0", api::PatchSide::T0 )
        .value( "T1", api::PatchSide::T1 );

    py::class_<topology::Cell>( m, "Cell" )
        .def_property_readonly(
            "dart",
            []( const topology::Cell& c ) { return c.dart().id(); },
            "Returns the id of a dart (or half edge) on the element." );

    py::class_<topology::Face, topology::Cell>( m, "Element" )
        .def( py::init( []( const topology::Dart::IndexType& id ) { return topology::Face( topology::Dart( id ) ); } ),
              "Constructs an element with the given dart id.",
              "dart"_a )
        .def( "__repr__", []( const topology::Face& c ) {
            std::stringstream ss;
            ss << c;
            return ss.str();
        } )
        .def( pybind11::self < pybind11::self )
        .def( pybind11::self == pybind11::self );

    py::class_<topology::Edge>( m, "Edge" )
        .def( py::init( []( const topology::Dart::IndexType& id ) { return topology::Edge( topology::Dart( id ) ); } ),
              "Constructs an edge with the given dart id.",
              "dart"_a )
        .def( "__repr__", []( const topology::Edge& c ) {
            std::stringstream ss;
            ss << c;
            return ss.str();
        } )
        .def( pybind11::self < pybind11::self )
        .def( pybind11::self == pybind11::self );

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
        .def(
            "evaluateBasisValuesAtParentPoint",
            &eval::SplineSpaceEvaluator::evaluateBasisValuesAtParentPoint,
            "Evaluates the extracted spline components at the currently localized parent (element-local) point. "
            "Vector-conforming and top-form components remain patch-parametric; no spatial Piola transformation is "
            "applied. The result has numFunctions(elem) rows and numVectorComponents columns." )
        .def(
            "evaluateBasis",
            &eval::SplineSpaceEvaluator::evaluateBasisValuesAtParentPoint,
            "Alias for evaluateBasisValuesAtParentPoint. The evaluation point is parent-local, while vector and "
            "top-form components remain patch-parametric." )
        .def(
            "evaluateBasisFirstDerivativesWrtParentCoordinates",
            &eval::SplineSpaceEvaluator::evaluateBasisFirstDerivativesWrtParentCoordinates,
            "Evaluates first derivatives of the extracted spline components with respect to parent coordinates. "
            "Columns are grouped by parent-coordinate direction, with one column per vector component in each group." )
        .def(
            "evaluateBasisFirstDerivativesWrtParametricCoordinates",
            &eval::SplineSpaceEvaluator::evaluateBasisFirstDerivativesWrtParametricCoordinates,
            "Evaluates first derivatives of the extracted spline components with respect to patch-parametric "
            "coordinates. Columns are grouped by parametric direction, with one column per vector component in each "
            "group." )
        .def(
            "evaluateFirstDerivatives",
            &eval::SplineSpaceEvaluator::evaluateBasisFirstDerivativesWrtParentCoordinates,
            "Alias for evaluateBasisFirstDerivativesWrtParentCoordinates. Derivatives are taken with respect to parent "
            "(element-local) coordinates." );

    const auto localize_parent_point = []( api::NavierStokesDiscretization& nsd, const Eigen::Vector2d& pt ) {
        const param::ParentPoint ppt( param::cubeDomain( 2 ), pt, { false, false, false, false } );
        nsd.getH1().localizeParentPoint( ppt );
        nsd.getHDIV().localizeParentPoint( ppt );
        nsd.getL2().localizeParentPoint( ppt );
        nsd.getGeometry().localizeParentPoint( ppt );
    };

    const auto parent_to_spatial_jacobian = []( const api::NavierStokesDiscretization& nsd ) {
        return nsd.getGeometry().evaluateParentToSpatialJacobian( nsd.geometryControlPoints() );
    };

    const auto parent_to_spatial_jacobian_determinant = []( const api::NavierStokesDiscretization& nsd ) {
        return nsd.getGeometry().evaluateParentToSpatialJacobian( nsd.geometryControlPoints() ).determinant();
    };

    const auto spatial_hdiv_basis = []( const api::NavierStokesDiscretization& nsd ) {
        return evaluateSpatialHDivBasisValues(
            nsd.getHDIV(),
            nsd.getGeometry(),
            nsd.geometryControlPoints() );
    };

    const auto spatial_hdiv_derivatives = []( const api::NavierStokesDiscretization& nsd ) {
        return evaluateSpatialHDivBasisFirstDerivatives(
            nsd.getHDIV(),
            nsd.getGeometry(),
            nsd.geometryControlPoints() );
    };

    const auto spatial_l2_basis = []( const api::NavierStokesDiscretization& nsd ) {
        return evaluateSpatialL2BasisValues(
            nsd.getL2(),
            nsd.getGeometry(),
            nsd.geometryControlPoints() );
    };

    const auto spatial_l2_derivatives = []( const api::NavierStokesDiscretization& nsd ) {
        return evaluateSpatialL2BasisFirstDerivatives(
            nsd.getL2(),
            nsd.getGeometry(),
            nsd.geometryControlPoints() );
    };

    py::class_<api::NavierStokesDiscretization>( m, "NavierStokesDiscretization" )
        .def_property_readonly(
            "H1", []( const api::NavierStokesDiscretization& d ) { return d.getH1(); }, "The H1 spline space." )
        .def_property_readonly(
            "HDIV",
            []( const api::NavierStokesDiscretization& d ) { return d.getHDIV(); },
            "The HDiv spline space, with two vector components." )
        .def_property_readonly(
            "L2", []( const api::NavierStokesDiscretization& d ) { return d.getL2(); }, "The L2 spline space." )
        .def_property_readonly(
            "geometry",
            []( const api::NavierStokesDiscretization& d ) { return d.getGeometry(); },
            "The spline-space evaluator used for geometry." )
        .def_property_readonly( "control_points",
                                &api::NavierStokesDiscretization::controlPoints,
                                "The Euclidean control points for the geometry." )
        .def_property_readonly( "geometry_control_points",
                                &api::NavierStokesDiscretization::geometryControlPoints,
                                "The coefficients used by the geometry evaluator. Rational geometries use homogeneous "
                                "control points." )
        .def_property_readonly( "is_rational",
                                &api::NavierStokesDiscretization::hasRationalGeometry,
                                "Whether the geometry is rational." )
        .def(
            "elements",
            []( const api::NavierStokesDiscretization& nsd ) {
                std::vector<topology::Face> out;
                out.reserve( cellCount( nsd.getH1().splineSpace().basisComplex().parametricAtlas().cmap(), 2 ) );
                iterateCellsWhile(
                    nsd.getH1().splineSpace().basisComplex().parametricAtlas().cmap(), 2, [&out]( const topology::Face& f ) {
                        out.push_back( f );
                        return true;
                    } );
                return out;
            },
            "A list of all the elements in the discretization." )
        .def(
            "localizeElement",
            []( api::NavierStokesDiscretization& nsd, const topology::Face& elem ) {
                nsd.getH1().localizeElement( elem );
                nsd.getHDIV().localizeElement( elem );
                nsd.getL2().localizeElement( elem );
                nsd.getGeometry().localizeElement( elem );
            },
            "Sets all future evaluations of the spline spaces to be on element elem, until localizeElement is called "
            "again.",
            "elem"_a )
        .def(
            "localizeParentPoint",
            localize_parent_point,
            "Localizes all evaluators at the given parent (element-local) coordinates in the currently localized "
            "element.",
            "pt"_a )
        .def(
            "localizePoint",
            localize_parent_point,
            "Alias for localizeParentPoint. The input is a parent (element-local) point, not a knot-space parametric "
            "point.",
            "pt"_a )
        .def(
            "mapping",
            []( const api::NavierStokesDiscretization& nsd ) {
                return nsd.getGeometry().evaluateManifold( nsd.geometryControlPoints() );
            },
            "Evaluates the geometric manifold map at the currently localized parent point and returns spatial "
            "coordinates." )
        .def(
            "parametricPoint",
            []( const api::NavierStokesDiscretization& nsd ) { return nsd.getGeometry().evaluateParametricPoint(); },
            "Returns the knot-space parametric coordinates corresponding to the currently localized parent point." )
        .def(
            "parentToSpatialJacobian",
            parent_to_spatial_jacobian,
            "Evaluates the parent-to-spatial geometry Jacobian dx/ds at the currently localized parent point." )
        .def(
            "jacobian",
            parent_to_spatial_jacobian,
            "Alias for parentToSpatialJacobian. Returns the parent-to-spatial geometry Jacobian dx/ds." )
        .def(
            "parentToSpatialJacobianDeterminant",
            parent_to_spatial_jacobian_determinant,
            "Evaluates det(dx/ds) for a square parent-to-spatial geometry Jacobian." )
        .def(
            "jacobianDeterminant",
            parent_to_spatial_jacobian_determinant,
            "Alias for parentToSpatialJacobianDeterminant. Returns det(dx/ds) for a square geometry Jacobian." )
        .def(
            "evaluateSpatialHDivBasisValues",
            spatial_hdiv_basis,
            "Evaluates the spatial H(div) basis using the contravariant Piola transformation from the patch-parametric "
            "frame. The basis is "
            "returned as a matrix with HDIV.numFunctions(elem) rows and two spatial-component columns, with each row "
            "corresponding to the same position in HDIV.connectivity(elem)." )
        .def(
            "piolaTransformedHDIVBasis",
            spatial_hdiv_basis,
            "Compatibility alias for evaluateSpatialHDivBasisValues." )
        .def(
            "evaluateSpatialHDivBasisFirstDerivatives",
            spatial_hdiv_derivatives,
            "Evaluates spatial derivatives of the H(div) basis. The result has HDIV.numFunctions(elem) rows and "
            "columns ordered as dv_x/dx, dv_y/dx, dv_x/dy, dv_y/dy." )
        .def(
            "piolaTransformedHDIVFirstDerivatives",
            spatial_hdiv_derivatives,
            "Compatibility alias for evaluateSpatialHDivBasisFirstDerivatives." )
        .def(
            "evaluateSpatialL2BasisValues",
            spatial_l2_basis,
            "Evaluates the spatial L2 top-form basis using its parametric-to-spatial Piola transformation. The result has "
            "L2.numFunctions(elem) rows and one column, with each row corresponding to the same position in "
            "L2.connectivity(elem)." )
        .def(
            "piolaTransformedL2",
            spatial_l2_basis,
            "Compatibility alias for evaluateSpatialL2BasisValues." )
        .def(
            "evaluateSpatialL2BasisFirstDerivatives",
            spatial_l2_derivatives,
            "Evaluates spatial derivatives of the L2 top-form basis. The result has L2.numFunctions(elem) rows and "
            "columns ordered as dp/dx, dp/dy." )
        .def(
            "piolaTransformedL2FirstDerivatives",
            spatial_l2_derivatives,
            "Compatibility alias for evaluateSpatialL2BasisFirstDerivatives." );

    py::class_<api::NavierStokesTPDiscretization, api::NavierStokesDiscretization>( m, "NavierStokesTPDiscretization" )
        .def( py::init<const basis::KnotVector&,
                       const basis::KnotVector&,
                       const size_t,
                       const size_t,
                       const Eigen::MatrixXd&>(),
              "Create a B-spline discretization for Navier Stokes problems, consisting of an H1 spline space with the "
              "given knot vectors, an HDiv spline space taking the H1 space as its primal basis, and an L2 spline "
              "space, which has reduced degree from H1 in both directions. The geometry is defined by the H1 space and "
              "the provided control points.",
              "knot_vec_s"_a,
              "knot_vec_t"_a,
              "degree_s"_a,
              "degree_t"_a,
              "control_points"_a )
        .def(
            "boundaryEdges",
            []( const api::NavierStokesTPDiscretization& nsd, const api::PatchSide& side ) {
                const topology::TPCombinatorialMap& tp_map = nsd.H1_ss.basisComplex().parametricAtlas().cmap();
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
            "boundaryEdges",
            []( const api::NavierStokesDiscretization& nsd ) {
                std::vector<topology::Edge> out;
                out.reserve( cellCount( nsd.cmapBdry(), 1 ) );
                iterateCellsWhile( nsd.cmapBdry(), 1, [&out]( const topology::Edge& e ) {
                    out.push_back( e );
                    return true;
                } );
                return out;
            },
            "A list of all edges on the boundary of the discretization" )
        .def(
            "boundaryPerpendicularHDivFuncs",
            []( const api::NavierStokesTPDiscretization& nsd, const api::PatchSide& side ) {
                const auto& component_bases = nsd.HDIV_ss.scalarTPBases();
                const bool is_S = side == api::PatchSide::S0 or side == api::PatchSide::S1;

                const util::IndexVec lengths = getTPLengths( is_S ? *component_bases.at( 0 ) : *component_bases.at( 1 ) );

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
        .def( "knotVectors", []( const api::NavierStokesTPDiscretization& nsd ) {
            const auto comps = tensorProductComponentSplines( nsd.H1_ss );
            return std::vector<basis::KnotVector>{ comps.at( 0 )->knotVector(), comps.at( 1 )->knotVector() };
        } );

    py::class_<api::NURBSNavierStokesTPDiscretization, api::NavierStokesTPDiscretization>(
        m, "NURBSNavierStokesTPDiscretization" )
        .def( py::init<const basis::KnotVector&,
                       const basis::KnotVector&,
                       const size_t,
                       const size_t,
                       const Eigen::MatrixXd&,
                       const Eigen::VectorXd&>(),
              "Create a tensor-product Navier-Stokes discretization with rational geometry. The analysis spaces remain "
              "polynomial; geometry is evaluated with homogeneous NURBS control points built from the Euclidean control "
              "points and weights.",
              "knot_vec_s"_a,
              "knot_vec_t"_a,
              "degree_s"_a,
              "degree_t"_a,
              "control_points"_a,
              "weights"_a )
        .def_property_readonly( "weights",
                                &api::NURBSNavierStokesTPDiscretization::weights,
                                "The rational geometry weights." );

    py::class_<api::NavierStokesHierarchicalDiscretization, api::NavierStokesDiscretization>( m, "NavierStokesHierarchicalDiscretization" )
        .def( py::init<const basis::KnotVector&,
                       const basis::KnotVector&,
                       const size_t,
                       const size_t,
                       const Eigen::MatrixXd&,
                       const std::vector<std::vector<std::pair<size_t, size_t>>>>(),
              "Create a B-spline discretization for Navier Stokes problems, consisting of an H1 spline space with the "
              "given knot vectors, an HDiv spline space taking the H1 space as its primal basis, and an L2 spline "
              "space, which has reduced degree from H1 in both directions. The geometry is defined by the H1 space and "
              "the provided control points.",
              "knot_vec_s"_a,
              "knot_vec_t"_a,
              "degree_s"_a,
              "degree_t"_a,
              "control_points"_a,
              "elems_to_refine"_a )
        .def(
            "boundaryEdges",
            []( const api::NavierStokesHierarchicalDiscretization& nsd, const api::PatchSide& side ) {
                const topology::HierarchicalTPCombinatorialMap& cmap = nsd.H1_ss.basisComplex().parametricAtlas().cmap();
                const topology::TPCombinatorialMap& base_cmap = *cmap.refinementLevels().front();
                std::vector<topology::Edge> out;

                const auto add_active_descendants = [&]( const topology::Edge& e ) {
                    const topology::Dart global_d = cmap.dartRanges().toGlobalDart( 0, e.dart() );
                    cmap.iterateLeafDescendants( global_d, [&]( const topology::Dart& leaf_d ) {
                        out.push_back( leaf_d );
                        return true;
                    } );
                };

                switch( side )
                {
                    case api::PatchSide::S0:
                    {
                        out.reserve( cellCount( base_cmap.lineCMap(), 1 ) );
                        const topology::Dart source_dart( 0 );
                        iterateDartsWhile( base_cmap.lineCMap(), [&]( const topology::Dart& line_dart ) {
                            add_active_descendants(
                                base_cmap.flatten( source_dart, line_dart, topology::TPCombinatorialMap::TPDartPos::DartPos3 ) );
                            return true;
                        } );
                        break;
                    }
                    case api::PatchSide::S1:
                    {
                        out.reserve( cellCount( base_cmap.lineCMap(), 1 ) );
                        const topology::Dart source_dart( base_cmap.sourceCMap().maxDartId() );
                        iterateDartsWhile( base_cmap.lineCMap(), [&]( const topology::Dart& line_dart ) {
                            add_active_descendants(
                                base_cmap.flatten( source_dart, line_dart, topology::TPCombinatorialMap::TPDartPos::DartPos1 ) );
                            return true;
                        } );
                        break;
                    }
                    case api::PatchSide::T0:
                    {
                        out.reserve( cellCount( base_cmap.sourceCMap(), 1 ) );
                        const topology::Dart line_dart( 0 );
                        iterateDartsWhile( base_cmap.sourceCMap(), [&]( const topology::Dart& source_dart ) {
                            add_active_descendants(
                                base_cmap.flatten( source_dart, line_dart, topology::TPCombinatorialMap::TPDartPos::DartPos0 ) );
                            return true;
                        } );
                        break;
                    }
                    case api::PatchSide::T1:
                    {
                        out.reserve( cellCount( base_cmap.sourceCMap(), 1 ) );
                        const topology::Dart line_dart( base_cmap.sourceCMap().maxDartId() );
                        iterateDartsWhile( base_cmap.sourceCMap(), [&]( const topology::Dart& source_dart ) {
                            add_active_descendants(
                                base_cmap.flatten( source_dart, line_dart, topology::TPCombinatorialMap::TPDartPos::DartPos2 ) );
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
            "boundaryEdges",
            []( const api::NavierStokesDiscretization& nsd ) {
                std::vector<topology::Edge> out;
                out.reserve( cellCount( nsd.cmapBdry(), 1 ) );
                iterateCellsWhile( nsd.cmapBdry(), 1, [&out]( const topology::Edge& e ) {
                    out.push_back( e );
                    return true;
                } );
                return out;
            },
            "A list of all edges on the boundary of the discretization" )
        .def(
            "boundaryPerpendicularHDivFuncs",
            []( const api::NavierStokesHierarchicalDiscretization& nsd, const api::PatchSide& side ) {
                const bool is_S = side == api::PatchSide::S0 or side == api::PatchSide::S1;

                const auto& scalar_hier_bases = nsd.HDIV_ss.scalarBases();
                const basis::HierarchicalTPSplineSpace& hier_comp = is_S ? *scalar_hier_bases.at( 0 ) : *scalar_hier_bases.at( 1 );

                const auto get_iter_dir = [&side]( const util::IndexVec& lengths ) -> SmallVector<std::variant<bool, size_t>, 3> {
                    switch( side )
                    {
                        case api::PatchSide::S0: return { size_t{ 0 }, true };
                        case api::PatchSide::S1: return { lengths.at( 0 ) - 1, true };
                        case api::PatchSide::T0: return { true, size_t{ 0 } };
                        case api::PatchSide::T1: return { true, lengths.at( 1 ) - 1 };
                    }
                };

                const size_t offset = is_S ? 0 : scalar_hier_bases.at( 0 )->numFunctions();
                const size_t num_levels = hier_comp.basisComplex().parametricAtlas().cmap().numLevels();

                std::vector<size_t> result; // LIKELY SLOW as I'm not reserving.
                for( size_t i = 0; i < num_levels; i++ )
                {
                    const basis::TPSplineSpace& comp = *hier_comp.refinementLevels().at( i );
                    const util::IndexVec lengths = getTPLengths( comp );
                    const auto iter_dir = get_iter_dir( lengths );
                    const auto& exop = hier_comp.levelExtractionOperators().at( i );

                    util::iterateTensorProduct( lengths, { 0, 1 }, iter_dir, [&]( const util::IndexVec& iv ) {
                        const size_t fid = util::flatten( iv, lengths );
                        for( Eigen::SparseMatrix<double>::InnerIterator it( exop, fid ); it; ++it )
                            result.push_back( it.row() + offset );
                    } );
                }

                std::ranges::sort( result );
                result.erase( std::ranges::unique( result ).end(), result.end() );
                return result;
            },
            "A list of all HDIV functions which are nonzero on and perpendicular to a given side of the discretization "
            "spline patch",
            "side"_a );

    m.def(
        "grevillePoints",
        []( const basis::KnotVector& kv_s, const basis::KnotVector& kv_t, const size_t degree_s, const size_t degree_t )
            -> Eigen::Matrix2Xd {
            return util::tensorProduct( { grevillePoints( kv_s, degree_s ), grevillePoints( kv_t, degree_t ) } ).transpose();
        },
        "Returns the greville points of the 2d B-spline patch defined by the given knot vectors and degrees.  These "
        "points, when used as control points, create a linear spatial spline geometry corresponding exactly to the "
        "parametic domain.",
        "kv_s"_a,
        "kv_t"_a,
        "degree_s"_a,
        "degree_t"_a );

    m.def(
        "globallyHRefine",
        []( const api::NavierStokesTPDiscretization& coarse_nsd,
            const size_t num_divisions,
            const double param_tol ) -> std::unique_ptr<api::NavierStokesTPDiscretization> {
            const auto component_bsplines = tensorProductComponentSplines( coarse_nsd.H1_ss );
            SmallVector<basis::KnotVector, 3> coarse_kvs;
            std::transform( component_bsplines.begin(),
                            component_bsplines.end(),
                            std::back_inserter( coarse_kvs ),
                            []( const auto& comp ) { return comp->knotVector(); } );

            SmallVector<basis::KnotVector, 3> fine_kvs;
            std::transform( coarse_kvs.begin(),
                            coarse_kvs.end(),
                            std::back_inserter( fine_kvs ),
                            [&num_divisions]( const auto& kv ) { return basis::nAdicRefine( kv, num_divisions ); } );

            util::IndexVec degrees;
            std::transform( component_bsplines.begin(),
                            component_bsplines.end(),
                            std::back_inserter( degrees ),
                            []( const auto& comp ) {
                                return comp->basisComplex().defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 );
                            } );

            const Eigen::MatrixXd fine_cpts =
                coarse_nsd.geometryControlPoints() * basis::refinementOp( coarse_kvs, fine_kvs, degrees, param_tol );

            if( coarse_nsd.hasRationalGeometry() )
                return std::make_unique<api::NURBSNavierStokesTPDiscretization>(
                    fine_kvs.at( 0 ), fine_kvs.at( 1 ), degrees.at( 0 ), degrees.at( 1 ), fine_cpts );

            return std::make_unique<api::NavierStokesTPDiscretization>(
                fine_kvs.at( 0 ), fine_kvs.at( 1 ), degrees.at( 0 ), degrees.at( 1 ), fine_cpts );
        },
        "Globally h-refines the given navier stokes discretization by evenly dividing each cell into num_divisions "
        "cells in each parametric dimension.",
        "coarse_nsd"_a,
        "num_divisions"_a,
        "parametric_tolerance"_a );

    m.def(
        "globallyPElevate",
        []( const api::NavierStokesTPDiscretization& coarse_nsd,
            const size_t degree_increase_s,
            const size_t degree_increase_t,
            const double param_tol ) -> std::unique_ptr<api::NavierStokesTPDiscretization> {
            const auto component_bsplines = basis::tensorProductComponentSplines( coarse_nsd.H1_ss );
            if( component_bsplines.size() != 2 )
                throw std::runtime_error( "globallyPElevate currently supports only tensor-product surface splines." );

            SmallVector<size_t, 3> target_degrees;
            target_degrees.push_back(
                component_bsplines.at( 0 )->basisComplex().defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 ) +
                degree_increase_s );
            target_degrees.push_back(
                component_bsplines.at( 1 )->basisComplex().defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 ) +
                degree_increase_t );

            const Eigen::SparseMatrix<double> elevation_op =
                basis::degreeElevationOp( coarse_nsd.H1_ss, target_degrees, param_tol );
            const Eigen::MatrixXd elevated_cpts = coarse_nsd.geometryControlPoints() * elevation_op;
            const SmallVector<basis::KnotVector, 3> elevated_kvs =
                basis::degreeElevatedKnotVectors( coarse_nsd.H1_ss, target_degrees );

            if( coarse_nsd.hasRationalGeometry() )
                return std::make_unique<api::NURBSNavierStokesTPDiscretization>(
                    elevated_kvs.at( 0 ), elevated_kvs.at( 1 ), target_degrees.at( 0 ), target_degrees.at( 1 ), elevated_cpts );

            return std::make_unique<api::NavierStokesTPDiscretization>(
                elevated_kvs.at( 0 ), elevated_kvs.at( 1 ), target_degrees.at( 0 ), target_degrees.at( 1 ), elevated_cpts );
        },
        "Globally degree-elevates the H1 tensor-product spline space of the given Navier-Stokes discretization.",
        "coarse_nsd"_a,
        "degree_increase_s"_a,
        "degree_increase_t"_a,
        "parametric_tolerance"_a );
}
