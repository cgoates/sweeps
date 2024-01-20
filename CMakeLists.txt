cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 23)
add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_LINKER_TYPE mold)
project(attempt)
add_executable( attempt main.cpp )

add_library( attemptlib
  Simplex.cpp
  VTKOutput.cpp
  AbaqusInput.cpp
  Logging.cpp
  Laplace.cpp
  SimplexUtilities.cpp
  Tracing.cpp
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
