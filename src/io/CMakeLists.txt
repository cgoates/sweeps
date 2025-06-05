add_library( io SHARED
  VTKOutput.cpp
  MeshInput.cpp
  QuadLayoutInput.cpp
  OBJOutput.cpp
)

target_link_libraries( io PUBLIC util topology basis )

target_include_directories( io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
