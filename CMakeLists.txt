cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 23)
add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_LINKER_TYPE mold)
add_compile_definitions("SRC_HOME=\"${CMAKE_SOURCE_DIR}\"")

project(sweep)

set(Catch2_DIR ${CMAKE_SOURCE_DIR}/../thirdparty/Catch2/install/lib/cmake/Catch2)

include_directories( ${CMAKE_SOURCE_DIR}/../thirdparty/eigen )

find_package(Eigen3 REQUIRED)

add_subdirectory(src)

find_package(Catch2 REQUIRED)

include(CTest)
include(Catch)
add_subdirectory(test)
