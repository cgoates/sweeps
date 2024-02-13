#pragma once
#include <string>
#include <SweepInput.hpp>
#include<map>

namespace io
{
    SweepInput loadINPFile( const std::string& filename,
                            const std::string& zero_bc_label,
                            const std::string& one_bc_label );

    void rewriteBaryCoordFile( const SimplicialComplex& mesh3d,
                            const std::string& obj_filename,
                            const std::string& bary_filename_in,
                            const std::string& bary_filename_out );

    std::vector<BarycentricPoint> loadBaryCoords( const std::string& bary_filename_in, const std::set<size_t>& max_edges );
}