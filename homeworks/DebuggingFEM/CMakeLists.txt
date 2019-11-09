
include(GoogleTest)

set(mastersolution_sources
  mastersolution/local_laplace_qfe.h
  mastersolution/local_laplace_qfe.cc
  mastersolution/main.cc
  mastersolution/qfe_interpolator.h
  mastersolution/qfe_interpolator.cc
  mastersolution/qfe_provider_tester.h
)

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

set(mastersolution_test_sources
  test/debugging_fem_test_mastersolution.cc
)

# Libraries

add_library(DebuggingFEM_mastersolution.solution ${mastersolution_sources})

# Executables

add_executable(DebuggingFEM_mastersolution ${mastersolution_sources})

add_executable(DebuggingFEM_test_mastersolution ${mastersolution_test_sources})


target_link_libraries(DebuggingFEM_mastersolution
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

target_link_libraries(DebuggingFEM_mastersolution.solution
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

target_link_libraries(DebuggingFEM_test_mastersolution 
  PUBLIC Eigen3::Eigen
  LF::lf.mesh.test_utils
  LF::lf.assemble
  GTest::gtest_main
  DebuggingFEM_mastersolution.solution
)

gtest_discover_tests(DebuggingFEM_test_mastersolution)
