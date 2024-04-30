add_library( topology SHARED
  TetMeshCombinatorialMap.cpp
  TriMeshCombinatorialMap.cpp
  CombinatorialMapMethods.cpp
  CombinatorialMapBoundary.cpp
  CombinatorialMapRestriction.cpp
  DelaunayTriangulation.cpp
  Dart.cpp
  Cell.cpp
  LevelSetCMap.cpp
)

target_link_libraries( topology PUBLIC util simplex )

target_include_directories( topology PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )