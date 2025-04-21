#pragma once
#include <vector>
#include <set>
#include <numbers>
#include <SimplicialComplex.hpp>

struct SweepInput
{
    SweepInput( const SimplicialComplex& m, const std::vector<bool>& z, const std::vector<bool>& o )
        : mesh( m ), zero_bcs( z ), one_bcs( o )
    {
        if( z.size() != m.points.size() or o.size() != m.points.size() )
            throw std::runtime_error( "Bad arguments to sweep input constructor" );
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