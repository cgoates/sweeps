add_subdirectory( util )
add_subdirectory( simplex )
add_subdirectory( topology )
add_subdirectory( param )
add_subdirectory( basis )
add_subdirectory( io )
add_subdirectory( reparam )

add_executable( sweep main.cpp )

target_link_libraries( sweep PUBLIC reparam )
target_include_directories( sweep PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )