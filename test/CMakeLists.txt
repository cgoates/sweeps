function(add_my_test arg)
    add_executable(${arg} ${arg}.cpp)
    target_link_libraries(${arg} PRIVATE Catch2::Catch2WithMain attemptlib)
    catch_discover_tests(${arg})
endfunction()

add_my_test(cgogn_tests)