add_subdirectory( util )
add_subdirectory( simplex )
add_subdirectory( topology )
add_subdirectory( io )
add_subdirectory( param )

add_executable( sweep main.cpp )

target_link_libraries( sweep PUBLIC param Eigen3::Eigen )
target_include_directories( sweep PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )