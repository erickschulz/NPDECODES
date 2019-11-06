include(GoogleTest)

set(mastersolution_sources
  mastersolution/handling_dofs.cc
  mastersolution/handling_dofs.h
  mastersolution/handling_dofs_main.cc
  )

set(test_mastersolution_sources 
  mastersolution/handling_dofs.cc
  mastersolution/handling_dofs.h
  test/handling_dofs_test_master.cc
  )

add_executable(HandlingDOFs_mastersolution ${mastersolution_sources})


add_executable(HandlingDOFs_test_mastersolution ${test_mastersolution_sources})

target_link_libraries(HandlingDOFs_mastersolution
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

target_link_libraries(HandlingDOFs_test_mastersolution
  GTest::gtest_main
  Eigen3::Eigen
  Boost::filesystem 
  Boost::system 
  Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  )

gtest_discover_tests(HandlingDOFs_test_mastersolution)
