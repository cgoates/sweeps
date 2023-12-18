cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 23)
add_compile_options(-Wall -Wextra -Wpedantic)
project(attempt)
add_executable( attempt main.cpp )

add_library( attemptlib
  SimplicialComplex.cpp
  TetMeshCombinatorialMap.cpp )

include_directories( ${CMAKE_SOURCE_DIR}/../eigen )

target_link_libraries( attempt PUBLIC attemptlib )
target_include_directories( attemptlib PUBLIC ${CMAKE_SOURCE_DIR} )
target_include_directories( attempt PUBLIC ${CMAKE_SOURCE_DIR} )
#include(ExternalProject)
#ExternalProject_Add( eigen SOURCE_DIR ${CMAKE_SOURCE_DIR}/../eigen)
#ExternalProject_Add( cgogn SOURCE_DIR ${CMAKE_SOURCE_DIR}/../CGOGN_3 )