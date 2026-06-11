#include <TraceExtraction.hpp>
#include <BasisComplex.hpp>
#include <CombinatorialMapMethods.hpp>
#include <IndexOperations.hpp>
#include <ParametricAtlas.hpp>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>

namespace basis
{
    namespace
    {
        util::IndexVec componentLengths( const SmallVector<size_t, 3>& degrees )
        {
            util::IndexVec lengths;
            std::transform( degrees.begin(), degrees.end(), std::back_inserter( lengths ), []( const size_t degree ) {
                return degree + 1;
            } );
            return lengths;
        }

        void checkSideAxis( const size_t dim, const ElementSide& side )
        {
            if( side.axis >= dim )
                throw std::invalid_argument( "ElementSide axis is incompatible with the parent basis dimension." );
        }
    }

    ElementSide elementSideFromId( const size_t side_id )
    {
        return ElementSide( side_id / 2, side_id % 2 == 1 );
    }

    SmallVector<size_t, 3> traceDegrees( const SmallVector<size_t, 3>& degrees, const ElementSide& side )
    {
        checkSideAxis( degrees.size(), side );

        SmallVector<size_t, 3> out;
        for( size_t i = 0; i < degrees.size(); i++ )
        {
            if( i != side.axis ) out.push_back( degrees.at( i ) );
        }
        return out;
    }

    std::vector<Eigen::Index> traceColumnIndices( const SmallVector<size_t, 3>& degrees,
                                                  const ElementSide& side )
    {
        checkSideAxis( degrees.size(), side );

        const util::IndexVec full_lengths = componentLengths( degrees );
        const util::IndexVec side_lengths = util::dropIndex( full_lengths, side.axis );
        std::vector<Eigen::Index> out;
        out.reserve( std::accumulate( side_lengths.begin(), side_lengths.end(), 1, std::multiplies<>() ) );

        util::iterateTensorProduct( side_lengths, [&]( const util::IndexVec& side_index ) {
            util::IndexVec full_index;
            for( size_t i = 0, side_i = 0; i < degrees.size(); i++ )
            {
                if( i == side.axis )
                    full_index.push_back( side.lower ? 0 : degrees.at( i ) );
                else
                    full_index.push_back( side_index.at( side_i++ ) );
            }
            out.push_back( static_cast<Eigen::Index>( util::flatten( full_index, full_lengths ) ) );
        } );

        return out;
    }

    std::vector<Eigen::Index> traceColumnIndices( const ParentBasis& pb, const ElementSide& side )
    {
        if( not param::isCartesian( pb.mParentDomain ) )
            throw std::runtime_error( "Trace extraction is only implemented for tensor-product parent bases." );
        checkSideAxis( param::dim( pb.mParentDomain ), side );

        if( numVectorComponents( pb ) == 1 ) return traceColumnIndices( degrees( pb ), side );

        std::vector<Eigen::Index> out;
        Eigen::Index component_offset = 0;
        for( const ParentBasis& component_basis : componentBases( pb ) )
        {
            const std::vector<Eigen::Index> component_cols = traceColumnIndices( component_basis, side );
            std::transform( component_cols.begin(),
                            component_cols.end(),
                            std::back_inserter( out ),
                            [&]( const Eigen::Index col ) { return component_offset + col; } );
            component_offset += static_cast<Eigen::Index>( numFunctions( component_basis ) );
        }

        return out;
    }

    Eigen::MatrixXd traceExtractionMatrix( const Eigen::MatrixXd& element_extraction,
                                           const ParentBasis& pb,
                                           const ElementSide& side )
    {
        if( element_extraction.cols() != static_cast<Eigen::Index>( numFunctions( pb ) ) )
            throw std::runtime_error( "Element extraction matrix column count does not match parent basis." );

        const std::vector<Eigen::Index> cols = traceColumnIndices( pb, side );
        Eigen::MatrixXd out( element_extraction.rows(), cols.size() );
        for( Eigen::Index j = 0; j < static_cast<Eigen::Index>( cols.size() ); j++ )
        {
            out.col( j ) = element_extraction.col( cols.at( j ) );
        }

        return out;
    }

    TraceSideData traceSideData( const SplineSpace& ss,
                                 const topology::Cell& element,
                                 const ElementSide& side,
                                 const double row_tol )
    {
        if( element.dim() != ss.basisComplex().parametricAtlas().cmap().dim() )
            throw std::runtime_error( "Trace side data must be requested from a top-dimensional element." );

        const ParentBasis pb = ss.basisComplex().parentBasis( element );
        const std::vector<FunctionId> element_connectivity = ss.connectivity( element );
        const Eigen::MatrixXd element_trace = traceExtractionMatrix( ss.extractionOperator( element ), pb, side );

        TraceSideData out{ element, side, {}, Eigen::MatrixXd(), {} };
        for( Eigen::Index row = 0; row < element_trace.rows(); row++ )
        {
            if( element_trace.row( row ).norm() > row_tol )
            {
                out.connectivity.push_back( element_connectivity.at( static_cast<size_t>( row ) ) );
                out.element_row_indices.push_back( static_cast<size_t>( row ) );
            }
        }

        out.extraction = Eigen::MatrixXd::Zero( out.connectivity.size(), element_trace.cols() );
        for( Eigen::Index trace_row = 0; trace_row < static_cast<Eigen::Index>( out.element_row_indices.size() );
             trace_row++ )
        {
            out.extraction.row( trace_row ) =
                element_trace.row( static_cast<Eigen::Index>( out.element_row_indices.at( trace_row ) ) );
        }

        return out;
    }

    std::optional<std::pair<topology::Cell, ElementSide>>
        adjacentElementSide( const topology::TPCombinatorialMap& cmap,
                             const topology::Cell& element,
                             const ElementSide& side )
    {
        if( element.dim() != cmap.dim() )
            throw std::runtime_error( "Adjacent element lookup requires a top-dimensional element." );
        checkSideAxis( cmap.dim(), side );

        const SmallVector<std::shared_ptr<const topology::CombinatorialMap1d>, 3> components =
            topology::tensorProductComponentCMaps( cmap );
        if( components.size() != cmap.dim() )
            throw std::runtime_error( "Adjacent element lookup requires a cube-like tensor-product map." );

        topology::FullyUnflattenedDart unflat = topology::unflattenFull( cmap, element.dart() );
        unflat.dart_pos =
            SmallVector<topology::TPCombinatorialMap::TPDartPos, 2>(
                cmap.dim() - 1, topology::TPCombinatorialMap::TPDartPos::DartPos0 );

        const size_t axis_length = topology::cellCount( *components.at( side.axis ), components.at( side.axis )->dim() );
        const size_t axis_index = unflat.unflat_darts.at( side.axis ).id();

        if( side.lower )
        {
            if( axis_index == 0 ) return std::nullopt;
            unflat.unflat_darts.at( side.axis ) = topology::Dart( axis_index - 1 );
        }
        else
        {
            if( axis_index + 1 >= axis_length ) return std::nullopt;
            unflat.unflat_darts.at( side.axis ) = topology::Dart( axis_index + 1 );
        }

        return std::pair<topology::Cell, ElementSide>{
            topology::Cell( topology::flattenFull( cmap, unflat ), cmap.dim() ), side.opposite() };
    }

    TraceInterfaceElement traceInterfaceElement( const SplineSpace& ss,
                                                 const topology::Cell& element,
                                                 const ElementSide& side,
                                                 const double row_tol )
    {
        const auto* tp_cmap =
            dynamic_cast<const topology::TPCombinatorialMap*>( &ss.basisComplex().parametricAtlas().cmap() );
        if( tp_cmap == nullptr )
            throw std::runtime_error( "traceInterfaceElement currently requires a tensor-product spline space." );

        TraceInterfaceElement out{ traceSideData( ss, element, side, row_tol ), std::nullopt };
        const std::optional<std::pair<topology::Cell, ElementSide>> neighbor =
            adjacentElementSide( *tp_cmap, element, side );
        if( neighbor.has_value() )
        {
            out.second = traceSideData( ss, neighbor->first, neighbor->second, row_tol );
        }

        return out;
    }

    TraceInterfaceElement traceInterfaceElement( const SplineSpace& first_ss,
                                                 const topology::Cell& first_element,
                                                 const ElementSide& first_side,
                                                 const SplineSpace& second_ss,
                                                 const topology::Cell& second_element,
                                                 const ElementSide& second_side,
                                                 const double row_tol )
    {
        return TraceInterfaceElement{
            traceSideData( first_ss, first_element, first_side, row_tol ),
            traceSideData( second_ss, second_element, second_side, row_tol ) };
    }
}
