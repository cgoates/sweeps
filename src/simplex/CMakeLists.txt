add_library( simplex SHARED
  Simplex.cpp
  SimplicialComplex.cpp
)

target_link_libraries( simplex PUBLIC util )

target_include_directories( simplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

