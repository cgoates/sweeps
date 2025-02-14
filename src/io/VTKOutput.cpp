#include <VTKOutput.hpp>
#include <fstream>
#include <Simplex.hpp>
#include <cassert>
#include <sstream>
#include <SplineSpace.hpp>
#include <BasisComplex.hpp>
#include <ParametricAtlas.hpp>
#include <CombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParentBasis.hpp>
#include <IndexOperations.hpp>
#include <SimplexUtilities.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace io
{
    inline std::string pieceHeader( const size_t n_points, const size_t n_cells )
    {
        std::ostringstream ss;
        ss <<R"STRING(<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">
  <UnstructuredGrid>
    <FieldData>
    </FieldData>
    <Piece NumberOfPoints=")STRING"
             << n_points << R"STRING(" NumberOfCells=")STRING" << n_cells << "\">";
        return ss.str();
    }

    inline std::string pieceFooter()
    {
        return "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>";
    }

    void VTKOutputObject::addVertexField( const std::string& label, const Eigen::Ref<const Eigen::MatrixXd> field )
    {
        mVertexFields.emplace( label, field );
    }

    void VTKOutputObject::addCellField( const std::string& label, const Eigen::Ref<const Eigen::MatrixXd> field )
    {
        mCellFields.emplace( label, field );
    }

    void outputSimplicialFieldToVTK( const VTKOutputObject& to_output, const std::string& filename )
    {
        const size_t n_simplices = to_output.complex().simplices.size();
        const size_t n_points = to_output.complex().points.size();

        std::ofstream file;
        file.open( filename );

        file << pieceHeader( n_points, n_simplices ) << "\n<CellData>";

        for( const auto& [label, data] : to_output.cellFields() )
        {
            file << R"STRING(
        <DataArray type="Float64" Name=")STRING" << label << R"STRING(" NumberOfComponents=")STRING"
                 << data.cols() << R"STRING(" format="ascii">)STRING" << std::endl;

            for( Eigen::Index row = 0; row < data.rows(); row++ )
            {
                for( Eigen::Index col = 0; col < data.cols(); col++ )
                {
                    file << data( row, col ) << " ";
                }
            }

            file << R"STRING(
        </DataArray>)STRING";
        }

        file << R"STRING(
      </CellData>
      <PointData>)STRING";

        for( const auto& [label, data] : to_output.vertexFields() )
        {
            file << R"STRING(
        <DataArray type="Float64" Name=")STRING" << label << R"STRING(" NumberOfComponents=")STRING"
                 << data.cols() << R"STRING(" format="ascii">)STRING" << std::endl;

            for( Eigen::Index row = 0; row < data.rows(); row++ )
            {
                for( Eigen::Index col = 0; col < data.cols(); col++ )
                {
                    file << data( row, col ) << " ";
                }
            }

            file << R"STRING(
        </DataArray>)STRING";
        }

        file << R"STRING(
      </PointData>
      <Points>
        <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">
)STRING";

        for( const Eigen::Vector3d& point : to_output.complex().points )
        {
            file << point( 0 ) << " " << point( 1 ) << " " << point( 2 ) << " ";
        }

        file << R"STRING(
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
)STRING";

        for( const Simplex& simplex : to_output.complex().simplices )
        {
            for( size_t vert_ii = 0; vert_ii <= simplex.dim(); vert_ii++ )
            {
                file << simplex.vertex( vert_ii ).id() << " ";
            }
        }

        file << R"STRING(
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
)STRING";

        size_t offset = 0;
        for( const Simplex& simplex : to_output.complex().simplices )
        {
            offset += simplex.dim() + 1;
            file << offset << " ";
        }
        file << R"STRING(
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
)STRING";

        for( const Simplex& simplex : to_output.complex().simplices )
        {
            switch( simplex.dim() )
            {
                // See documentation at https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf page 9
                case 0: file << 1 << " "; break;
                case 1: file << 3 << " "; break;
                case 2: file << 5 << " "; break;
                case 3: file << 10 << " "; break;
                default: break;
            }
        }

        file << R"STRING(
        </DataArray>
      </Cells>
)STRING" << pieceFooter();

        file.close();
    }

    SmallVector<size_t, 3> degreesOfParentBasis( const basis::ParentBasis& pb )
    {
        SmallVector<size_t, 3> out;
        for( const auto& g : pb.mBasisGroups )
        {
            for( const auto& p : g.degrees )
                out.push_back( p );
        }
        return out;
    }

    void outputBezierMeshToVTK( const basis::SplineSpace& ss,
                                const MatrixX3dMax& geom,
                                const std::string& filename )
    {
        outputBezierMeshToVTK( BezierOutputObject( ss, geom ), filename );
    }

    void outputBezierMeshToVTK( const BezierOutputObject& bo, const std::string& filename )
    {
        outputPartialBezierMeshToVTK(
            bo, filename, [&]( const std::function<void( const topology::Cell& )>& callback ) {
                iterateCellsWhile( bo.ss().basisComplex().parametricAtlas().cmap(),
                                   bo.ss().basisComplex().parametricAtlas().cmap().dim(),
                                   [&]( const topology::Cell& c ) {
                                       callback( c );
                                       return true;
                                   } );
            } );
    }

    void outputPartialBezierMeshToVTK( const basis::SplineSpace& ss,
                                       const MatrixX3dMax& geom,
                                       const std::string& filename,
                                       const std::function<void( const std::function<void( const topology::Cell& )>& )>& cell_iterator )
    {
        outputPartialBezierMeshToVTK( BezierOutputObject( ss, geom ), filename, cell_iterator );
    }

    Eigen::MatrixXd degreeElevationMatrix( const size_t source_degree, const size_t target_degree )
    {
        if( target_degree < source_degree ) throw std::invalid_argument( "Cannot degree elevate to a lower degree" );

        Eigen::MatrixXd out = Eigen::MatrixXd::Identity( source_degree + 1, source_degree + 1 );
        for( size_t p = source_degree + 1; p <= target_degree; p++ )
        {
            Eigen::MatrixXd elev_p = Eigen::MatrixXd::Zero( p, p + 1 );
            for( size_t n = 0; n < p; n++ )
            {
                elev_p( n, n ) = double( p - n ) / double( p );
                elev_p( n, n + 1 ) = double( n + 1 ) / double( p );
            }
            out = out * elev_p;
        }

        return out;
    }

    Eigen::MatrixXd degreeElevate( const Eigen::MatrixXd& in, const basis::ParentBasis& source_basis, const basis::ParentBasis& target_basis )
    {
        if( source_basis == target_basis ) return in;
        if( numVectorComponents( target_basis ) != 1 ) throw std::invalid_argument( "Cannot output to vector basis" );

        if( numVectorComponents( source_basis ) > 1 )
        {
            if( in.cols() > 1 ) throw std::invalid_argument( "Vector valued bases cannot be expanded with vector valued control points" );
            const SmallVector<basis::ParentBasis, 3> component_pbs = componentBases( source_basis );

            Eigen::MatrixXd out( in.rows(), component_pbs.size() );
            for( size_t i = 0; i < component_pbs.size(); i++ )
            {
                out.col( i ) = degreeElevate( in, component_pbs.at( i ), target_basis );
            }

            return out;
        }
        else
        {
            if( not isCartesian( source_basis.mParentDomain ) or not isCartesian( target_basis.mParentDomain ) )
                throw std::invalid_argument( "Only cartesian bases are supported for degree elevation" );

            const SmallVector<size_t, 3> source_degrees = degreesOfParentBasis( source_basis );
            const SmallVector<size_t, 3> target_degrees = degreesOfParentBasis( target_basis );

            Eigen::MatrixXd elevation_matrix = Eigen::MatrixXd::Ones( 1, 1 );

            for( size_t i = 0; i < source_degrees.size(); i++ )
            {
                elevation_matrix = Eigen::kroneckerProduct( degreeElevationMatrix( source_degrees.at( i ), target_degrees.at( i ) ), elevation_matrix ).eval();
            }

            return elevation_matrix.transpose() * in;
        }
    }

    void outputPartialBezierMeshToVTK( const BezierOutputObject& bo,
                                       const std::string& filename,
                                       const std::function<void( const std::function<void( const topology::Cell& )>& )>& cell_iterator )
    {
        const auto& ss = bo.ss();
        const auto& geom = bo.geom();
        const topology::CombinatorialMap& cmap = ss.basisComplex().parametricAtlas().cmap();
        const size_t param_dim = cmap.dim();
        const size_t spatial_dim = geom.cols();

        std::ostringstream degrees_string;
        std::ostringstream points_string;
        std::ostringstream connectivity_string;
        std::ostringstream offsets_string;
        std::ostringstream types_string;
        std::map<std::string, std::ostringstream> field_strings;
        for( const auto& [label, field] : bo.bezierFields() )
        {
            field_strings.emplace( label, std::ostringstream() );
        }

        size_t n_points = 0;
        size_t n_cells = 0;
        cell_iterator( [&]( const topology::Cell& c ) {
            n_cells++;
            const basis::ParentBasis pb = ss.basisComplex().parentBasis( c );
            SmallVector<size_t, 3> degrees = degreesOfParentBasis( pb );
            SmallVector<size_t, 3> tp_lengths;
            for( const size_t&  p : degrees ) tp_lengths.push_back( p + 1 );
            while( degrees.size() < 3 ) degrees.push_back( 0 );
            degrees_string << degrees.at( 0 ) << " " << degrees.at( 1 ) << " " << degrees.at( 2 ) << std::endl;

            const Eigen::MatrixXd ex_op = ss.extractionOperator( c );
            const std::vector<basis::FunctionId> connect = ss.connectivity( c );
            const Eigen::MatrixXd bez_points = ex_op.transpose() * geom( connect, Eigen::all );
            
            util::iterateVTKTPOrdering( tp_lengths, [&] ( const size_t row_ii ) {
                if( spatial_dim == 3 )
                    points_string << bez_points( row_ii, 0 ) << " " << bez_points( row_ii, 1 ) << " " << bez_points( row_ii, 2 ) << std::endl;
                else if( spatial_dim == 2 )
                    points_string << bez_points( row_ii, 0 ) << " " << bez_points( row_ii, 1 ) << " 0" << std::endl;
                else if( spatial_dim == 1 )
                    points_string << bez_points( row_ii, 0 ) << " 0 0" << std::endl;
            } );

            for( Eigen::Index i = 0; i < bez_points.rows(); i++ )
                connectivity_string << ( n_points + i ) << " ";
            n_points += bez_points.rows();
            offsets_string << n_points << " ";

            switch( param_dim )
            {
                // See documentation at
                // https://www.kitware.com/main/wp-content/uploads/2020/03/Implementation-of-rational-Be%CC%81zier-cells-into-VTK-Report.pdf
                // page 3
                case 0: types_string << 1 << " "; break;
                case 1: types_string << 75 << " "; break;
                case 2: types_string << 77 << " "; break;
                case 3: types_string << 79 << " "; break;
                default: break;
            }

            // Add the fields
            for( const auto& [label, field] : bo.bezierFields() )
            {
                const basis::SplineSpace& field_ss = field.ss.has_value() ? field.ss.value().get() : ss;
                const Eigen::MatrixXd field_ex_op = field.ss.has_value() ? field_ss.extractionOperator( c ) : ex_op;
                const std::vector<basis::FunctionId> field_conn = field.ss.has_value() ? field_ss.connectivity( c ) : connect;
                const Eigen::MatrixXd field_vals = field_ex_op.transpose() * field.geom( field_conn, Eigen::all );

                const Eigen::MatrixXd field_bez_points = degreeElevate( field_vals, field_ss.basisComplex().parentBasis( c ), pb );

                util::iterateVTKTPOrdering( tp_lengths, [&] ( const size_t row_ii ) {
                    for( size_t i = 0; i < field_bez_points.cols() - 1; i++ )
                        field_strings.at( label ) << field_bez_points( row_ii, i ) << " ";
                    field_strings.at( label ) << field_bez_points( row_ii, field_bez_points.cols() - 1 ) << std::endl;
                } );

            }
        } );

        std::ofstream file;
        file.open( filename );

        file << pieceHeader( n_points, n_cells ) << "\n<CellData HigherOrderDegrees=\"HigherOrderDegrees\">";

        file << R"STRING(
        <DataArray type="UInt8" Name="HigherOrderDegrees" NumberOfComponents="3" format="ascii">)STRING" << std::endl;

        file << degrees_string.str();

        file << R"STRING(
        </DataArray>
      </CellData>
      <PointData>
)STRING";

        for( const auto& [label, field] : field_strings )
        {
            file << R"STRING(<DataArray type="Float64" Name=")STRING" << label << R"STRING(" NumberOfComponents=")STRING"
                 << bo.bezierFields().at( label ).geom.cols() << R"STRING(" format="ascii">)STRING" << std::endl;
            file << field.str();
            file << "</DataArray>";
        }

        file << R"STRING(
      </PointData>
      <Points>
        <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">
)STRING";

        file << points_string.str();

        file << R"STRING(
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
)STRING";

        file << connectivity_string.str();

        file << R"STRING(
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
)STRING";

        file << offsets_string.str();

        file << R"STRING(
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
)STRING";

        file << types_string.str();

        file << R"STRING(
        </DataArray>
      </Cells>
)STRING" << pieceFooter();

        file.close();
    }

    void outputEdges( const topology::CombinatorialMap& cmap,
                      const VertexPositionsFunc& positions,
                      const std::vector<topology::Edge>& edges,
                      const std::string& filename )
    {
        SimplicialComplex path;

        if( not edges.empty() )
        {
            path.points.push_back( positions( topology::Vertex( edges.front().dart() ) ) );
            for( const auto& edge : edges )
            {
                const size_t temp = path.points.size();
                path.points.push_back( positions( topology::Vertex( phi( cmap, 1, edge.dart() ).value() ) ) );
                path.simplices.emplace_back( temp - 1, temp );
            }
        }

        io::VTKOutputObject output( path );
        io::outputSimplicialFieldToVTK( output, filename );
    }

    void outputCMap( const topology::CombinatorialMap& cmap,
                     const VertexPositionsFunc& positions,
                     const std::string& filename )
    {
        SimplicialComplex out;
        std::map<topology::Dart::IndexType, size_t> vert_ids;
        iterateCellsWhile( cmap, 0, [&]( const auto& vert ) {
            vert_ids.emplace( lowestDartId( cmap, vert ), out.points.size() );
            out.points.push_back( positions( vert ) );
            return true;
        } );

        iterateCellsWhile( cmap, 2, [&]( const auto& f ) {
            const topology::Dart& d = f.dart();
            const VertexId v1 = vert_ids.at( lowestDartId( cmap, topology::Vertex( d ) ) );
            const VertexId v2 = vert_ids.at( lowestDartId( cmap, topology::Vertex( phi( cmap, 1, d ).value() ) ) );
            const VertexId v3 = vert_ids.at( lowestDartId( cmap, topology::Vertex( phi( cmap, -1, d ).value() ) ) );
            out.simplices.emplace_back( v1, v2, v3 );
            return true;
        } );

        io::VTKOutputObject output( out );
        io::outputSimplicialFieldToVTK( output, filename );
    }

    void outputDualFace( const topology::CombinatorialMap& cmap,
                         const VertexPositionsFunc& positions,
                         const topology::Edge& e,
                         const std::string& postfix )
    {
        SimplicialComplex dual;
        SimplicialComplex tets;
        SimplicialComplex edges;

        const Eigen::Vector3d edge_mid = 0.5 * ( positions( topology::Vertex( e.dart() ) ) +
                                                 positions( topology::Vertex( phi( cmap, 1, e.dart() ).value() ) ) );

        topology::Dart curr_d = e.dart();
        if( boundaryAdjacent( cmap, e ) )
        {
            // iterate phi {2,3} until there is no phi({2,3})
            std::optional<topology::Dart> maybe_next = std::nullopt;
            do
            {
                maybe_next = phi( cmap, {2,3}, curr_d );
                if( maybe_next.has_value() )
                    curr_d = maybe_next.value();
            } while( maybe_next.has_value() );
        }

        const topology::Dart start_dart = curr_d;

        Eigen::Vector3d last_circumcenter = circumcenter( tetOfVolume( cmap, positions, topology::Volume( curr_d ) ) );
        do
        {
            const topology::Volume tet( curr_d );
            addTetNoDuplicateChecking( tets, cmap, positions, tet );
            addEdgeNoDuplicateChecking( edges, cmap, positions, topology::Edge( curr_d ) );

            const std::optional<topology::Dart> maybe_next_d = phi( cmap, {3,2}, curr_d );
            if( maybe_next_d.has_value() )
            {
                const topology::Dart next_d = maybe_next_d.value();
                const Eigen::Vector3d this_circumcenter = circumcenter( tetOfVolume( cmap, positions, topology::Volume( next_d ) ) );
                dual.points.push_back( edge_mid );
                dual.points.push_back( last_circumcenter );
                dual.points.push_back( this_circumcenter );
                dual.simplices.emplace_back( dual.points.size() - 3, dual.points.size() - 2, dual.points.size() - 1 );
                last_circumcenter = this_circumcenter;
                curr_d = next_d;
            }
            else
            {
                break;
            }
        } while ( curr_d != start_dart );

        io::VTKOutputObject output( dual );
        io::outputSimplicialFieldToVTK( output, "dual" + postfix + ".vtu" );

        io::VTKOutputObject output_tets( tets );
        io::outputSimplicialFieldToVTK( output_tets, "tetsneardual" + postfix + ".vtu" );

        io::VTKOutputObject output_edges( edges );
        io::outputSimplicialFieldToVTK( output_edges, "edgesneardual" + postfix + ".vtu" );
    }
} // namespace io