add_executable( navier_infsup NavierStokesInfSup.cpp )

target_link_libraries( navier_infsup PUBLIC spline_api io eval basis util )
target_include_directories( navier_infsup PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )