cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 23)
add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_LINKER_TYPE mold)
add_compile_definitions("SRC_HOME=\"${CMAKE_SOURCE_DIR}\"")

set(Catch2_DIR ${CMAKE_SOURCE_DIR}/../thirdparty/Catch2/install/lib/cmake/Catch2)


project(sweep)
add_executable( sweep main.cpp )

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

include_directories( ${CMAKE_SOURCE_DIR}/../thirdparty/eigen )

find_package(Eigen3 REQUIRED)

target_link_libraries( sweeplib PUBLIC Eigen3::Eigen )

target_link_libraries( sweep PUBLIC sweeplib Eigen3::Eigen )
target_include_directories( sweeplib PUBLIC ${CMAKE_SOURCE_DIR} )
target_include_directories( sweep PUBLIC ${CMAKE_SOURCE_DIR} )

find_package(Catch2 REQUIRED)

include(CTest)
include(Catch)
add_subdirectory(test)
