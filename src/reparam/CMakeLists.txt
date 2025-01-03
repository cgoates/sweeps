add_library( reparam SHARED
  Laplace.cpp
  Tracing.cpp
  Foliation.cpp
  CriticalPoints.cpp
)

target_link_libraries( reparam PUBLIC util topology mapping io )

target_include_directories( reparam PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )