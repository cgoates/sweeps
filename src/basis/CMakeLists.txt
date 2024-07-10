add_library( basis SHARED
  ParentBasis.cpp
  BasisComplex.cpp
  BasisComplex1d.cpp
  TPBasisComplex.cpp
  DivConfBasisComplex.cpp
  SplineSpace.cpp
  GenericSplineSpace.cpp
  TPSplineSpace.cpp
  KnotVector.cpp
)

target_link_libraries( basis PUBLIC util topology param )

target_include_directories( basis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )