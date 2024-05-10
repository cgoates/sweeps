add_library( reparam SHARED
  Laplace.cpp
  Tracing.cpp
  SimplexUtilities.cpp
)

target_link_libraries( reparam PUBLIC util topology io )

target_include_directories( reparam PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )