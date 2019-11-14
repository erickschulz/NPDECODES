
include(GoogleTest)


set(mastersolution_sources
  mastersolution/lin_fe_react_diff.cc
  mastersolution/lin_fe_react_diff.h
  mastersolution/main.cc
)

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

set(test_sources
  test/lin_fe_react_diff_tests.cc
)


add_library(LinFeReactDiff.mastersolution ${mastersolution_sources})


add_executable(LinFeReactDiff_mastersolution ${mastersolution_sources})
add_executable(LinFeReactDiff_test_mastersolution ${test_sources})


target_link_libraries(LinFeReactDiff_mastersolution
  PUBLIC Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)


target_link_libraries(LinFeReactDiff.mastersolution
  PUBLIC Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  )

target_link_libraries(LinFeReactDiff_test_mastersolution GTest::gtest_main LinFeReactDiff.mastersolution)

gtest_discover_tests(LinFeReactDiff_test_mastersolution)

