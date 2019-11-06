include(GoogleTest)


set(mastersolution_sources
mastersolution/gradprojection.h
mastersolution/main.cc
)
set(master_test_sources
  mastersolution/test.cc
)

add_executable(ProjectionOntoGradients_mastersolution ${mastersolution_sources})

target_link_libraries(ProjectionOntoGradients_mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

#Unit tests
add_executable(ProjectionOntoGradients_test_mastersolution ${master_test_sources})

add_library(ProjectionOntoGradients_master.solution ${mastersolution_sources})

target_link_libraries(ProjectionOntoGradients_master.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  )



target_link_libraries(ProjectionOntoGradients_test_mastersolution
GTest::gtest_main 
ProjectionOntoGradients_master.solution
)


gtest_discover_tests(ProjectionOntoGradients_test_mastersolution)



