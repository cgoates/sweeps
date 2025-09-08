find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
pkg_check_modules(PETSC REQUIRED PETSc)
pkg_check_modules(SLEPC REQUIRED slepc)

add_executable( navier_infsup NavierStokesInfSup.cpp SLEPcUtils.cpp )

target_link_libraries( navier_infsup PUBLIC spline_api io eval basis util MPI::MPI_CXX ${SLEPC_LIBRARIES} )
target_include_directories( navier_infsup PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MPI_CXX_INCLUDE_PATH}
  ${PETSC_INCLUDE_DIRS}
  ${SLEPC_INCLUDE_DIRS}
)

target_link_directories( navier_infsup PRIVATE
    ${SLEPC_LIBRARY_DIRS}
)


add_executable( maxwell_eigenproblem MaxwellEigenproblem.cpp SLEPcUtils.cpp DeRhamComplexDiscretization.cpp )

target_link_libraries( maxwell_eigenproblem PUBLIC spline_api io eval basis util MPI::MPI_CXX ${SLEPC_LIBRARIES} )
target_include_directories( maxwell_eigenproblem PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MPI_CXX_INCLUDE_PATH}
  ${PETSC_INCLUDE_DIRS}
  ${SLEPC_INCLUDE_DIRS}
)

target_link_directories( maxwell_eigenproblem PRIVATE
    ${SLEPC_LIBRARY_DIRS}
)