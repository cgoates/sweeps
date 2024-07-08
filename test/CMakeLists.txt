function(add_my_test arg)
    add_executable(${arg} ${arg}.cpp)
    target_link_libraries(${arg} PRIVATE Catch2::Catch2WithMain reparam basis )

    target_include_directories( ${arg} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
    catch_discover_tests(${arg})
endfunction()

Include(FetchContent)

set( FETCHCONTENT_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../deps")

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.0
)

FetchContent_MakeAvailable(Catch2)

include(CTest)
include(Catch)

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
add_my_test(Foliation_test)
add_my_test(TPCombinatorialMap_test)
add_my_test(TPParametricAtlas_test)
add_my_test(SplineSpace1d_test)