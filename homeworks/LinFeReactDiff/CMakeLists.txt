
include(GoogleTest)

set(template_sources
  templates/lin_fe_react_diff.cc
  templates/lin_fe_react_diff.h
  templates/main.cc
)

set(mysolution_sources
  mysolution/lin_fe_react_diff.cc
  mysolution/lin_fe_react_diff.h
  mysolution/main.cc
)

set(mastersolution_sources
  mastersolution/lin_fe_react_diff.cc
  mastersolution/lin_fe_react_diff.h
  mastersolution/main.cc
)

set(test_sources
  test/lin_fe_react_diff_tests.cc
)

add_library(LinFeReactDiff.solution ${mysolution_sources})

add_library(LinFeReactDiff.mastersolution ${mastersolution_sources})

add_executable(LinFeReactDiff_template ${template_sources})

add_executable(LinFeReactDiff_mysolution ${mysolution_sources})

add_executable(LinFeReactDiff_mastersolution ${mastersolution_sources})

add_executable(LinFeReactDiff_test_mysolution ${test_sources})

add_executable(LinFeReactDiff_test_mastersolution ${test_sources})


target_link_libraries(LinFeReactDiff_template
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(LinFeReactDiff_mysolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(LinFeReactDiff_mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(LinFeReactDiff.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
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
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  )

target_link_libraries(LinFeReactDiff_test_mysolution GTest::main LinFeReactDiff.solution)
target_link_libraries(LinFeReactDiff_test_mastersolution GTest::main LinFeReactDiff.mastersolution)


gtest_discover_tests(LinFeReactDiff_test_mysolution)
gtest_discover_tests(LinFeReactDiff_test_mastersolution)

