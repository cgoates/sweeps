add_library( util SHARED
  Logging.cpp
  SmallQueue.cpp
  SmallVector.cpp
)

add_dependencies( util eigen )

target_include_directories( util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )