#include <MFEMOutput.hpp>
#include <fstream>
#include <MultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <GlobalCellMarker.hpp>
#include <SmallQueue.hpp>

using namespace io;
using namespace basis;
using namespace topology;

constexpr std::array<topology::Dart::IndexType, 4> verts_2d = { 0, 1, 2, 3 };
constexpr std::array<topology::Dart::IndexType, 8> verts_3d = { 0, 6, 12, 18, 23, 5, 11, 17 };

std::vector<std::string> patchVertexAdjacencies( const MultiPatchCombinatorialMap& block_layout, const IndexingFunc& vert_ids )
{
    using namespace topology;
    std::vector<std::string> out;
    out.reserve( block_layout.constituents().size() );

    const auto iterate_vert_dart_ids = [&]( const std::function<void(const Dart::IndexType)>& callback ) {
        if( block_layout.dim() == 2 )
        {
            for( const auto& v : verts_2d )
            {
                callback( v );
            }
        }
        else
        {
            for( const auto& v : verts_3d )
            {
                callback( v );
            }
        }
    };
    for( size_t i = 0; i < block_layout.constituents().size(); i++ )
    {
        std::ostringstream oss;
        bool first = true;
        iterate_vert_dart_ids( [&]( const Dart::IndexType& did ) {
            const Vertex vert( did );
            const size_t vid = vert_ids( Vertex( block_layout.toGlobalDart( i, vert.dart() ) ) );
            if( first )
            {
                first = false;
                oss << vid;
            }
            else
            {
                oss << " " << vid;
            }
        } );
        out.push_back( oss.str() );
    }

    return out;
}

bool floodEquivalentEdges( const MultiPatchCombinatorialMap& block_layout, const Edge& start_edge,
                          const std::function<bool(const Edge&,const bool)>& callback )
{
    const uint dim = block_layout.dim();
    GlobalDartMarker m( block_layout );
    GrowableQueue<std::pair<Dart, bool>, 300> dart_queue;
    dart_queue.push( { start_edge.dart(), false } );
    m.mark( start_edge.dart() );

    const auto add_to_queue = [&m,&dart_queue]( const std::optional<Dart>& d, const bool reverse = false ) {
        if( d.has_value() and not m.isMarked( d.value() ) )
        {
            dart_queue.push( { d.value(), reverse } );
            m.mark( d.value() );
        }
    };
    while( not dart_queue.empty() )
    {
        const auto [curr_d, reverse] = dart_queue.pop();
        if( not callback( curr_d, reverse ) ) return false;
        if( dim == 3 )
        {
            add_to_queue( phi( block_layout, { 1, 1, 2 }, curr_d ) );
            add_to_queue( phi( block_layout, { 3, 2 }, curr_d ) );
        }
        else
        {
            const auto next_d = phi( block_layout, {1,1,2}, curr_d );
            if( not next_d.has_value() ) add_to_queue( phi( block_layout, {1, 1}, curr_d ), true );
            else add_to_queue( next_d );
        }
    }
    return true;
}

namespace io
{
void outputMultiPatchSplinesToMFEM( const basis::MultiPatchSplineSpace& mpss,
                                    const Eigen::Ref<const Eigen::MatrixXd> geom,
                                    const std::string& filename )
{
    outputMultiPatchSplinesToMFEM(
        mpss,
        geom,
        []( const size_t ) { return 1; }, // patch_to_block_id
        []( const MultiPatchCombinatorialMap& block_layout, const topology::Cell& c ) {
            const TPCombinatorialMap::TPDartPos pos = std::get<2>(
                block_layout.constituents().front()->unflatten( block_layout.toLocalDart( c.dart() ).second ) );
            const size_t bdry_group = pos == TPCombinatorialMap::TPDartPos::DartPos0   ? 1
                                      : pos == TPCombinatorialMap::TPDartPos::DartPos5 ? 2
                                                                                       : 3;
            return bdry_group;
        }, // bdry_cell_to_block_id
        filename );
}

void outputMultiPatchSplinesToMFEM( const basis::MultiPatchSplineSpace& mpss,
                                    const Eigen::Ref<const Eigen::MatrixXd> geom,
                                    const std::function<size_t( const size_t )>& patch_to_block_id,
                                    const std::function<size_t( const topology::MultiPatchCombinatorialMap&,
                                                                const topology::Cell& )>& bdry_cell_to_block_id,
                                    const std::string& filename )
{
    std::ofstream file;
    file.open( filename );

    const uint dim = mpss.basisComplex().parametricAtlas().cmap().dim();
    const MultiPatchCombinatorialMap block_layout = blockLayout( mpss.basisComplex().parametricAtlas().cmap() );

    const auto block_layout_vert_ids = [&]() {
        std::map<Dart::IndexType, size_t> ids;
        size_t id = 0;
        block_layout.iterateCellsWhile( 0, [&]( const Vertex& v ) {
            ids.emplace( lowestDartId( block_layout, v ), id++ );
            return true;
        } );
        return [&block_layout,ids]( const Vertex& v ) -> size_t {
            return ids.at( lowestDartId( block_layout, v ) );
        };
    }();

    file << "MFEM NURBS mesh v1.0\n\n";
    file << "dimension\n" << dim << "\n\n";

    file << "elements\n" << mpss.subSpaces().size() << "\n";
    const auto patch_adjacencies = patchVertexAdjacencies( block_layout, block_layout_vert_ids );
    const int elem_type = ( dim == 2 ) ? 3 : 5;
    const int bdry_type = ( dim == 2 ) ? 1 : 3;
    for( size_t patch_ii = 0; patch_ii < mpss.subSpaces().size(); patch_ii++ )
    {
        const size_t block_id = patch_to_block_id( patch_ii );
        file << block_id << " " << elem_type << " " << patch_adjacencies.at( patch_ii ) << "\n";
    }

    file << R"STRING(

boundary
)STRING";
    const CombinatorialMapBoundary bdry( block_layout );
    file << cellCount( bdry, dim - 1 ) << "\n";

    iterateCellsWhile( block_layout, dim - 1, [&]( const Cell& c ) {
        if( not boundaryAdjacent( block_layout, c ) )
        {
            return true; // skip non-boundary faces
        }
        const size_t bdry_group = bdry_cell_to_block_id( block_layout, c );

        file << bdry_group << " " << bdry_type;
        iterateDartsOfRestrictedCell( block_layout, c, dim - 1, [&]( const Dart& d ) {
            const size_t vert_id = block_layout_vert_ids( Vertex( d ) );
            file << " " << vert_id;
            return true;
        } );
        file << "\n";

        return true;
    } );

    file << "\nedges\n" << cellCount( block_layout, 1 ) << "\n";

    const auto add_edge = [&]( const size_t i, const Edge& edge, const bool reverse ) {
        const size_t vert_id1 = block_layout_vert_ids( Vertex( edge.dart() ) );
        const size_t vert_id2 = block_layout_vert_ids( Vertex( topology::phi( block_layout, 1, edge.dart() ).value() ) );
        if( reverse ) file << i << " " << vert_id2 << " " << vert_id1 << "\n";
        else file << i << " " << vert_id1 << " " << vert_id2 << "\n";
    };

    const SmallVector<Dart, 3> param_dim_darts = ( dim == 2 ) ? SmallVector<Dart, 3>{ Dart( 0 ), Dart( 3 ) }
                                                        : SmallVector<Dart, 3>{ Dart( 0 ), Dart( 19 ), Dart( 2 ) };

    GlobalCellMarker marker( block_layout, 1 );
    size_t kv_id = 0;
    for( size_t i = 0; i < mpss.subSpaces().size(); ++i )
    {
        for( const Dart patch_d : param_dim_darts )
        {
            const Edge edge( block_layout.toGlobalDart( i, patch_d ) );
            if( marker.isMarked( edge ) ) continue;

            floodEquivalentEdges( block_layout, edge, [&]( const Edge& e, const bool reverse ) {
                if( marker.isMarked( e ) ) return true; // already marked
                marker.mark( block_layout, e );
                add_edge( kv_id, e, reverse );
                return true;
            } );

            kv_id++;
        }
    }

    file << "\nvertices\n";

    file << cellCount( block_layout, 0 ) << "\n\npatches\n\n";

    const auto& func_ids = mpss.functionIdMap();
    for( size_t i = 0; i < mpss.subSpaces().size(); ++i )
    {
        const auto& patch = mpss.subSpaces().at( i );
        file << "knotvectors\n" << dim << "\n";
        const auto patch_components = tensorProductComponentSplines( *patch );
        for( const auto& comp : patch_components )
        {
            file << comp->basisComplex().defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 ) << " "
                 << comp->numFunctions();
            const auto& kv = comp->knotVector();
            for( size_t i = 0; i < kv.size(); ++i )
            {
                file << " " << kv.knot( i );
            }
            file << "\n";
        }
        file << "\ndimension\n" << dim << "\n\n";

        file << "controlpoints\n";
        for( const auto& func_id : func_ids.at( i ) )
        {
            constexpr int default_weight = 1;
            file << geom.col( func_id.id() ).transpose() << " " << default_weight << "\n";
        }
        file << "\n\n";
    }

    file.close();
}


}
