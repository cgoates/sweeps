function(add_my_test arg)
    add_executable(${arg} ${arg}.cpp)
    target_link_libraries(${arg} PRIVATE Catch2::Catch2WithMain reparam basis eval )

    target_include_directories( ${arg} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      catch_discover_tests(${arg})
    else()
      catch_discover_tests(${arg} TEST_SPEC "~[no_ubuntu]")
    endif()
endfunction()

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.0
)

FetchContent_MakeAvailable(Catch2)

include(CTest)
include(Catch)

add_my_test(CustomEigen_test)
add_my_test(tet_map_tests)
add_my_test(tri_map_tests)
add_my_test(boundary_map_tests)
add_my_test(levelset_map_tests)
add_my_test(simplex_utility_tests)
add_my_test(laplace_tests)
add_my_test(tracing_tests)
add_my_test(macaroni_vertex_tracing_test)
add_my_test(macaroni_cell_tracing_test)
add_my_test(macaroni_interior_face_tracing_test)
add_my_test(TriangleMeshCircleMapping_test)
add_my_test(TriangleMeshMapping_test)
add_my_test(Foliation_test)
add_my_test(TPCombinatorialMap_test)
add_my_test(MultiPatchCombinatorialMap_test)
add_my_test(TPParametricAtlas_test)
add_my_test(SplineSpace1d_test)
add_my_test(MultiPatchSplineSpace_test)
add_my_test(DivConfTPSplineSpace_test)
add_my_test(SplineSpaceEvaluator_test)
add_my_test(CoonsPatch_test)
target_link_libraries( CoonsPatch_test PUBLIC fitting)
add_my_test(IndexOperations_test)
add_my_test(HierarchicalTPCMap_test)
add_my_test(HierarchicalTPPAtlas_test)
add_my_test(HierarchicalTPSplineSpace_test)
add_my_test(HierarchicalMultiPatchCMap_test)
add_my_test(HierarchicalMultiPatchSplineSpace_test)
add_my_test(FitToFoliation_test)
target_link_libraries( FitToFoliation_test PUBLIC fitting)
add_my_test(Dijkstra_test)
add_my_test(CutCombinatorialMap_test)
add_my_test(NonDiskFoliation_test)
target_link_libraries(NonDiskFoliation_test PUBLIC fitting)
add_my_test(CriticalPoints_test)
add_my_test(SphereFoliation_test)
target_link_libraries(SphereFoliation_test PUBLIC fitting)
add_my_test(MultiPatchDecomposition_test)

add_test(
    NAME python_import_test
    COMMAND python3 -c "import sweeps; import splines; print('Import successful')"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/src/api
)