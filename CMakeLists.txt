cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 23)
add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_LINKER_TYPE mold)
add_compile_definitions("SRC_HOME=\"${CMAKE_SOURCE_DIR}\"")

set(Meshb_DIR ${CMAKE_SOURCE_DIR}/../thirdparty/cgogn_build/lib/cmake/Meshb)
set(cgogn_core_DIR ${CMAKE_SOURCE_DIR}/../thirdparty/cgogn_build/lib/cmake/cgogn_core)
set(cgogn_geometry_DIR ${CMAKE_SOURCE_DIR}/../thirdparty/cgogn_build/lib/cmake/cgogn_geometry)
set(cgogn_io_DIR ${CMAKE_SOURCE_DIR}/../thirdparty/cgogn_build/lib/cmake/cgogn_io)
set(Catch2_DIR ${CMAKE_SOURCE_DIR}/../thirdparty/Catch2/install/lib/cmake/Catch2)


project(attempt)
add_executable( attempt main.cpp )

add_library( attemptlib
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
)

include_directories( ${CMAKE_SOURCE_DIR}/../thirdparty/eigen )
include_directories( ${CMAKE_SOURCE_DIR}/../thirdparty/CGoGN_3 )

find_package(Eigen3 REQUIRED)
find_package(cgogn_core REQUIRED)
find_package(cgogn_io REQUIRED)
find_dependency(Meshb REQUIRED)

target_link_libraries( attemptlib PUBLIC cgogn::io cgogn::core )

target_link_libraries( attempt PUBLIC attemptlib Eigen3::Eigen cgogn::io cgogn::core )
target_include_directories( attemptlib PUBLIC ${CMAKE_SOURCE_DIR} )
target_include_directories( attempt PUBLIC ${CMAKE_SOURCE_DIR} )

find_package(Catch2 REQUIRED)

include(CTest)
include(Catch)
add_subdirectory(test)
