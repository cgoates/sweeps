add_library( io SHARED
  VTKOutput.cpp
  AbaqusInput.cpp
  QuadLayoutInput.cpp
)

target_link_libraries( io PUBLIC util topology )

target_include_directories( io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )