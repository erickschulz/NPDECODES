include(GoogleTest)

set(mastersolution_sources
  mastersolution/handling_dofs.cc
  mastersolution/handling_dofs.h
  mastersolution/handling_dofs_main.cc
  )

set(templates_sources
  templates/handling_dofs.cc
  templates/handling_dofs.h
  templates/handling_dofs_main.cc
  )

set(mysolution_sources
  mysolution/handling_dofs.cc
  mysolution/handling_dofs.h
  mysolution/handling_dofs_main.cc
  )

set(test_mysolution_sources 
  mysolution/handling_dofs.cc
  mysolution/handling_dofs.h
  test/handling_dofs_test.cc
  )

set(test_mastersolution_sources 
  mastersolution/handling_dofs.cc
  mastersolution/handling_dofs.h
  test/handling_dofs_test_master.cc
  )

add_executable(HandlingDOFs_mastersolution ${mastersolution_sources})

add_executable(HandlingDOFs_templates ${templates_sources})

add_executable(HandlingDOFs_mysolution ${mysolution_sources})

add_executable(HandlingDOFs_test_mysolution ${test_mysolution_sources})

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

target_link_libraries(HandlingDOFs_templates
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

target_link_libraries(HandlingDOFs_mysolution
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

target_link_libraries(HandlingDOFs_test_mysolution
  GTest::main
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

target_link_libraries(HandlingDOFs_test_mastersolution
  GTest::main
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

gtest_discover_tests(HandlingDOFs_test_mysolution)
gtest_discover_tests(HandlingDOFs_test_mastersolution)
