add_library( topology SHARED
  TetMeshCombinatorialMap.cpp
  TriMeshCombinatorialMap.cpp
  CombinatorialMapMethods.cpp
  CombinatorialMapBoundary.cpp
  CombinatorialMapRestriction.cpp
  DelaunayTriangulation.cpp
  CutCombinatorialMap.cpp
  Dart.cpp
  DartRange.cpp
  Cell.cpp
  LevelSetCMap.cpp
  ReversedCombinatorialMap.cpp
  CombinatorialMap1d.cpp
  TPCombinatorialMap.cpp
  MultiPatchCombinatorialMap.cpp
  HierarchicalTPCombinatorialMap.cpp
  HierarchicalMultiPatchCombinatorialMap.cpp
  SimplexUtilities.cpp
  Dijkstra.cpp
)

target_link_libraries( topology PUBLIC util simplex )

target_include_directories( topology PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )