add_library( fitting SHARED
  CoonsPatch.cpp
  LeastSquaresFitting.cpp
)

target_link_libraries( fitting PUBLIC util topology param basis eval )

target_include_directories( fitting PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )