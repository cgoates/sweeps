#pragma once
#include <vector>
#include <set>
#include <numbers>
#include <Eigen/Dense>
#include <SimplicialComplex.hpp>

struct SweepInput
{
    SweepInput( const SimplicialComplex& m, const std::vector<bool>& z, const std::vector<bool>& o )
        : mesh( m ), zero_bcs( z ), one_bcs( o )
    {
        if( z.size() != m.points.size() or o.size() != m.points.size() )
            throw( "Bad arguments to sweep input constructor" );
    }

    static SweepInput fromSets( const SimplicialComplex& m, const std::set<VertexId>& z, const std::set<VertexId>& o )
    {
        std::vector<bool> zero_bcs( m.points.size(), false );
        std::vector<bool> one_bcs( m.points.size(), false );
        for( const VertexId& vid : z ) zero_bcs.at( vid.id() ) = true;
        for( const VertexId& vid : o ) one_bcs.at( vid.id() ) = true;
        return SweepInput( m, zero_bcs, one_bcs );
    }

    SimplicialComplex mesh;
    std::vector<bool> zero_bcs;
    std::vector<bool> one_bcs;
};

class SweepInputTestCases
{
    public:
    static SweepInput twoTets()
    {
        std::vector<Simplex> simplices;
        simplices.push_back( { 0, 1, 2, 3 } );
        simplices.push_back( { 1, 4, 2, 3 } );

        std::vector<Eigen::Vector3d> points;
        points.push_back( { 0, 0, 0 } );
        points.push_back( { 1, 0, 0 } );
        points.push_back( { 0, 1, 0 } );
        points.push_back( { 0.5, 0.5, 0.5 } );
        points.push_back( { 1, 1, 0 } );

        return SweepInput::fromSets( { simplices, points }, { 0 }, { 1 } );
    }

    static SweepInput twelveTetCube()
    {
        std::vector<Simplex> simplices;
        simplices.push_back( { 0, 1, 2, 8 } );
        simplices.push_back( { 1, 3, 2, 8 } );
        simplices.push_back( { 0, 5, 1, 8 } );
        simplices.push_back( { 0, 4, 5, 8 } );
        simplices.push_back( { 2, 3, 7, 8 } );
        simplices.push_back( { 6, 2, 7, 8 } );
        simplices.push_back( { 1, 5, 3, 8 } );
        simplices.push_back( { 5, 7, 3, 8 } );
        simplices.push_back( { 2, 6, 4, 8 } );
        simplices.push_back( { 0, 2, 4, 8 } );
        simplices.push_back( { 5, 4, 7, 8 } );
        simplices.push_back( { 7, 4, 6, 8 } );

        std::vector<Eigen::Vector3d> points;
        points.push_back( { 0, 0, 0 } );
        points.push_back( { 1, 0, 0 } );
        points.push_back( { 0, 1, 0 } );
        points.push_back( { 1, 1, 0 } );
        points.push_back( { 0, 0, 1 } );
        points.push_back( { 1, 0, 1 } );
        points.push_back( { 0, 1, 1 } );
        points.push_back( { 1, 1, 1 } );
        points.push_back( { 0.5, 0.5, 0.5 } );

        return SweepInput::fromSets( { simplices, points }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } );
    }

    static SweepInput refinedCube( const Eigen::Vector3i& n_refinements )
    {
        // 5 tet cube

        std::vector<Simplex> simplices_template_1;
        simplices_template_1.push_back( { 0, 1, 2, 8 } );
        simplices_template_1.push_back( { 1, 3, 2, 8 } );
        simplices_template_1.push_back( { 0, 5, 1, 8 } );
        simplices_template_1.push_back( { 0, 4, 5, 8 } );
        simplices_template_1.push_back( { 2, 3, 7, 8 } );
        simplices_template_1.push_back( { 6, 2, 7, 8 } );
        simplices_template_1.push_back( { 1, 5, 3, 8 } );
        simplices_template_1.push_back( { 5, 7, 3, 8 } );
        simplices_template_1.push_back( { 2, 6, 4, 8 } );
        simplices_template_1.push_back( { 0, 2, 4, 8 } );
        simplices_template_1.push_back( { 5, 4, 6, 8 } );
        simplices_template_1.push_back( { 7, 5, 6, 8 } );

        const Eigen::Vector3d spacing = n_refinements.cast<double>().cwiseInverse();
        std::vector<Eigen::Vector3d> points;
        for( int k = 0; k <= n_refinements( 2 ); k++ )
        {
            for( int j = 0; j <= n_refinements( 1 ); j++ )
            {
                for( int i = 0; i <= n_refinements( 0 ); i++ )
                {
                    points.emplace_back( i * spacing( 0 ), j * spacing( 1 ), k * spacing( 2 ) );
                }
            }
        }
        const int n_cube_corners = points.size();

        const Eigen::Vector3d half_spacing = spacing * 0.5;
        for( int k = 0; k < n_refinements( 2 ); k++ )
        {
            for( int j = 0; j < n_refinements( 1 ); j++ )
            {
                for( int i = 0; i < n_refinements( 0 ); i++ )
                {
                    points.emplace_back( ( 1 + 2 * i ) * half_spacing( 0 ),
                                         ( 1 + 2 * j ) * half_spacing( 1 ),
                                         ( 1 + 2 * k ) * half_spacing( 2 ) );
                }
            }
        }

        const auto flatten_vert = [&n_refinements]( const Eigen::Vector3i& vert ) {
            return vert( 0 ) + ( n_refinements( 0 ) + 1 ) * vert( 1 ) +
                   ( n_refinements( 0 ) + 1 ) * ( n_refinements( 1 ) + 1 ) * vert( 2 );
        };

        const auto flatten_hex = [&n_refinements]( const Eigen::Vector3i& hex ) {
            return hex( 0 ) + n_refinements( 0 ) * hex( 1 ) + n_refinements( 0 ) * n_refinements( 1 ) * hex( 2 );
        };

        std::vector<Simplex> simplices;
        for( int k = 0; k < n_refinements( 2 ); k++ )
        {
            for( int j = 0; j < n_refinements( 1 ); j++ )
            {
                for( int i = 0; i < n_refinements( 0 ); i++ )
                {
                    // Find the points in the grid that correspond to the ijkth cube
                    const std::array<int, 9> local_verts = { flatten_vert( { i, j, k } ),
                                                             flatten_vert( { i + 1, j, k } ),
                                                             flatten_vert( { i, j + 1, k } ),
                                                             flatten_vert( { i + 1, j + 1, k } ),
                                                             flatten_vert( { i, j, k + 1 } ),
                                                             flatten_vert( { i + 1, j, k + 1 } ),
                                                             flatten_vert( { i, j + 1, k + 1 } ),
                                                             flatten_vert( { i + 1, j + 1, k + 1 } ),
                                                             n_cube_corners + flatten_hex( { i, j, k } ) };

                    for( const Simplex& simp : ( ( k % 2 == 0 ) ? simplices_template_1 : simplices_template_1 ) )
                        simplices.push_back( { VertexId( local_verts.at( simp.vertex( 0 ).id() ) ),
                                               VertexId( local_verts.at( simp.vertex( 1 ).id() ) ),
                                               VertexId( local_verts.at( simp.vertex( 2 ).id() ) ),
                                               VertexId( local_verts.at( simp.vertex( 3 ).id() ) ) } );
                }
            }
        }

        std::set<VertexId> zero_bcs;
        std::set<VertexId> one_bcs;
        for( int j = 0; j <= n_refinements( 1 ); j++ )
        {
            for( int i = 0; i <= n_refinements( 0 ); i++ )
            {
                zero_bcs.insert( flatten_vert( { i, j, 0 } ) );
                one_bcs.insert( flatten_vert( { i, j, n_refinements( 2 ) } ) );
            }
        }

        return SweepInput::fromSets( { simplices, points }, zero_bcs, one_bcs );
    }

    static SweepInput bentRefinedCube( const Eigen::Vector3i& n_refinements )
    {
        using namespace std::numbers;
        const SweepInput cube = refinedCube( n_refinements );

        std::vector<Eigen::Vector3d> new_points;
        std::transform(
            cube.mesh.points.begin(),
            cube.mesh.points.end(),
            std::back_inserter( new_points ),
            [&]( const Eigen::Vector3d& pt ) {
                return Eigen::Vector3d( 0.1 * pt( 0 ),
                                        ( 1 + 0.1 * pt( 1 ) ) * sin( pi * pt( 2 ) / 2 ),
                                        ( 1 + 0.1 * pt( 1 ) ) * cos( pi * pt( 2 ) / 2 ) );
            } );

        const SweepInput out( { cube.mesh.simplices, new_points }, cube.zero_bcs, cube.one_bcs );

        return out;
    }
};