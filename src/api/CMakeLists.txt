Include(FetchContent)

set( FETCHCONTENT_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../deps" )

find_package( Python COMPONENTS Interpreter Development REQUIRED )
include_directories( ${PYTHON_INCLUDE_DIRS} )

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11
  GIT_TAG        v2.13.1
)

FetchContent_MakeAvailable( pybind11 )

pybind11_add_module( splines SplinesAPI.cpp NavierStokesDiscretization.cpp )
target_link_libraries( splines PUBLIC basis eval )
target_include_directories( splines PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

pybind11_add_module( sweeps SweepAPI.cpp SweepAPIMethods.cpp )
target_link_libraries( sweeps PUBLIC basis eval io reparam )
target_include_directories( sweeps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

pybind11_add_module( eigenfunction EigenfunctionAPI.cpp )
target_include_directories( eigenfunction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )