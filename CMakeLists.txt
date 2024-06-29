cmake_minimum_required(VERSION 3.24)
set(CMAKE_CXX_STANDARD 23)
add_compile_options(-Wall -Wextra -Wpedantic)
add_compile_definitions("SRC_HOME=\"${CMAKE_SOURCE_DIR}\"")

project(sweep)

include(ExternalProject)

set(EIGEN_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen-install")

ExternalProject_Add(
    eigen
    URL  https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/deps"
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen-src"
    BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen-build"
    INSTALL_DIR "${EIGEN_INSTALL_DIR}"
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=Release
)

file( MAKE_DIRECTORY ${EIGEN_INSTALL_DIR}/include/eigen3 )

include_directories( SYSTEM ${EIGEN_INSTALL_DIR}/include/eigen3 )

add_subdirectory(src)

enable_testing()
add_subdirectory(test)
