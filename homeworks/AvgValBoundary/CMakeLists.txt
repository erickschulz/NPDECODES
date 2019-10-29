
include(GoogleTest)

set(template_sources
  templates/avg_val_boundary.cc
  templates/comp_gal_mat.h
  templates/main.cc
)

set(mysolution_sources
  mysolution/avg_val_boundary.cc
  mysolution/comp_gal_mat.h
  mysolution/main.cc
)

set(mastersolution_sources
  mastersolution/avg_val_boundary.cc
  mastersolution/comp_gal_mat.h
  mastersolution/main.cc
)

set(test_sources
  test/avg_val_boundary_tests.cc
)

set(master_test_sources
  test/avg_val_boundary_tests_master.cc
)

add_library(AvgValBoundary.solution ${mysolution_sources})

add_library(AvgValBoundary.mastersolution ${mastersolution_sources})

add_executable(AvgValBoundary_template ${template_sources})

add_executable(AvgValBoundary_mysolution ${mysolution_sources})

add_executable(AvgValBoundary_mastersolution ${mastersolution_sources})

add_executable(AvgValBoundary_test_mysolution ${test_sources})

add_executable(AvgValBoundary_test_mastersolution ${master_test_sources})


target_link_libraries(AvgValBoundary_template
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

target_link_libraries(AvgValBoundary_mysolution
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

target_link_libraries(AvgValBoundary_mastersolution
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

target_link_libraries(AvgValBoundary.solution
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

target_link_libraries(AvgValBoundary.mastersolution
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

target_link_libraries(AvgValBoundary_test_mysolution GTest::main AvgValBoundary.solution)
target_link_libraries(AvgValBoundary_test_mastersolution GTest::main AvgValBoundary.mastersolution)


gtest_discover_tests(AvgValBoundary_test_mysolution)
gtest_discover_tests(AvgValBoundary_test_mastersolution)
