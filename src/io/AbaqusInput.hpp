#pragma once
#include <string>
#include <SweepInput.hpp>

namespace io
{
    SweepInput loadINPFile( const std::string& filename,
                            const std::string& zero_bc_label,
                            const std::string& one_bc_label );
}