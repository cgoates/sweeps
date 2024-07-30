add_library( eval SHARED
  ParentBasisEval.cpp
  SplineSpaceEvaluator.cpp
)

target_link_libraries( eval PUBLIC util topology param basis )

target_include_directories( eval PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )