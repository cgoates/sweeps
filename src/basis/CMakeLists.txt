add_library( basis SHARED
  ParentBasis.cpp
  BasisComplex.cpp
  BasisComplex1d.cpp
  TPBasisComplex.cpp
  MultiPatchBasisComplex.cpp
  DivConfBasisComplex.cpp
  SplineSpace.cpp
  BSplineSpace1d.cpp
  TPSplineSpace.cpp
  MultiPatchSplineSpace.cpp
  DivConfTPSplineSpace.cpp
  KnotVector.cpp
  HierarchicalTPBasisComplex.cpp
)

target_link_libraries( basis PUBLIC util topology param )

target_include_directories( basis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )