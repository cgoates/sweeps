add_library( util
  Logging.cpp
  SmallQueue.cpp
  SmallVector.cpp
)

target_link_libraries( util PUBLIC Eigen3::Eigen )

target_include_directories( util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )