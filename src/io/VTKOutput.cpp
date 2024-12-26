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
        const topology::CombinatorialMap& cmap = ss.basisComplex().parametricAtlas().cmap();
        const size_t param_dim = cmap.dim();
        const size_t spatial_dim = geom.cols();
        const size_t n_cells = cellCount( cmap, param_dim );

        std::ostringstream degrees_string;
        std::ostringstream points_string;
        std::ostringstream connectivity_string;
        std::ostringstream offsets_string;
        size_t n_points = 0;
        iterateCellsWhile( cmap, param_dim, [&]( const topology::Cell& c ) {
            SmallVector<size_t, 3> degrees = degreesOfParentBasis( ss.basisComplex().parentBasis( c ) );
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

            return true;
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

        iterateCellsWhile( cmap, param_dim, [&]( const auto& ) {
            switch( param_dim )
            {
                // See documentation at
                // https://www.kitware.com/main/wp-content/uploads/2020/03/Implementation-of-rational-Be%CC%81zier-cells-into-VTK-Report.pdf
                // page 3
                case 0: file << 1 << " "; break;
                case 1: file << 75 << " "; break;
                case 2: file << 77 << " "; break;
                case 3: file << 79 << " "; break;
                default: break;
            }
            return true;
        } );

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
        const auto lowest_dart_id = [&cmap]( const topology::Vertex& v ) {
            topology::Dart::IndexType min_d = v.dart().id();
            iterateDartsOfCell( cmap, v, [&]( const topology::Dart& d ) {
                min_d = std::min( d.id(), min_d );
                return true;
            } );
            return min_d;
        };
        SimplicialComplex out;
        std::map<topology::Dart::IndexType, size_t> vert_ids;
        iterateCellsWhile( cmap, 0, [&]( const auto& vert ) {
            vert_ids.emplace( lowest_dart_id( vert ), out.points.size() );
            out.points.push_back( positions( vert ) );
            return true;
        } );

        iterateCellsWhile( cmap, 2, [&]( const auto& f ) {
            const topology::Dart& d = f.dart();
            const VertexId v1 = vert_ids.at( lowest_dart_id( d ) );
            const VertexId v2 = vert_ids.at( lowest_dart_id( phi( cmap, 1, d ).value() ) );
            const VertexId v3 = vert_ids.at( lowest_dart_id( phi( cmap, -1, d ).value() ) );
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