add_library( mapping SHARED
  GeometricMapping.cpp
  TriangleMeshCircleMapping.cpp
  SimplexUtilities.cpp
)

target_link_libraries( mapping PUBLIC util topology param )

target_include_directories( mapping PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )