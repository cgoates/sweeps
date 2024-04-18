function(add_my_test arg)
    add_executable(${arg} ${arg}.cpp)
    target_link_libraries(${arg} PRIVATE Catch2::Catch2WithMain sweeplib)

    target_include_directories( ${arg} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
    catch_discover_tests(${arg})
endfunction()

add_my_test(tet_map_tests)
add_my_test(tri_map_tests)
add_my_test(boundary_map_tests)
add_my_test(simplex_utility_tests)
add_my_test(laplace_tests)
add_my_test(tracing_tests)