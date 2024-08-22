add_library( param SHARED
  ParentDomain.cpp
  ParentPoint.cpp
  ParametricAtlas.cpp
  ParametricAtlas1d.cpp
  TriangleParametricAtlas.cpp
  TPParametricAtlas.cpp
  MultiPatchParametricAtlas.cpp
)

target_link_libraries( param PUBLIC util topology )

target_include_directories( param PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )