add_library( param SHARED
  Laplace.cpp
  Tracing.cpp
  SimplexUtilities.cpp
)

target_link_libraries( param PUBLIC util topology io )

target_include_directories( param PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )