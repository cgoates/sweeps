add_library( basis SHARED
  ParentBasis.cpp
  BasisComplex.cpp
  SplineSpace.cpp
)

target_link_libraries( basis PUBLIC util topology param )

target_include_directories( basis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )