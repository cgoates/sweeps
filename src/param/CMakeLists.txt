add_library( param SHARED
  ParentDomain.cpp
  ParametricAtlas.cpp
  TriangleParametricAtlas.cpp
)

target_link_libraries( param PUBLIC util topology )

target_include_directories( param PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )