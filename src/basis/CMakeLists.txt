add_library( basis SHARED
  LocalBasis.cpp
  BasisMesh.cpp
  SplineSpace.cpp
)

target_link_libraries( basis PUBLIC util topology param )

target_include_directories( basis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )