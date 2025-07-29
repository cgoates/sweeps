add_library( fitting SHARED
  CoonsPatch.cpp
  LeastSquaresFitting.cpp
  HexMeshOptimization.cpp
)

FetchContent_Declare(
  HexOpt
  GIT_REPOSITORY https://github.com/cgoates/HexOpt.git
  GIT_TAG main
)

FetchContent_MakeAvailable(HexOpt)

target_link_libraries( fitting PUBLIC util topology param basis eval HexOpt )

target_include_directories( fitting PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )