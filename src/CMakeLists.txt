add_library( sweeplib
  Simplex.cpp
  VTKOutput.cpp
  AbaqusInput.cpp
  QuadLayoutInput.cpp
  Logging.cpp
  Laplace.cpp
  SimplexUtilities.cpp
  Tracing.cpp
  TetMeshCombinatorialMap.cpp
  CombinatorialMapMethods.cpp
  CombinatorialMapBoundary.cpp
  CombinatorialMapRestriction.cpp
  Dart.cpp
  Cell.cpp
)

target_link_libraries( sweeplib PUBLIC Eigen3::Eigen )

target_include_directories( sweeplib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )


add_executable( sweep main.cpp )

target_link_libraries( sweep PUBLIC sweeplib Eigen3::Eigen )
target_include_directories( sweep PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )