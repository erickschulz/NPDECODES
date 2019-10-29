include(GoogleTest)

set(template_sources
  templates/test_problem.h
  templates/test_problem.cc
  templates/utilities.cc
)

set(mysolution_sources
  mysolution/test_problem.h
  mysolution/test_problem.cc
  mysolution/utilities.cc

)

set(mastersolution_sources
  mastersolution/test_problem.h
  mastersolution/test_problem.cc
  mastersolution/utilities.cc
)

set(test_sources
  test/test_problem_tests.cc
)

add_library(TestProblem.solution ${mysolution_sources})

add_library(TestProblem.mastersolution ${mastersolution_sources})

add_executable(TestProblem_template ${template_sources})

add_executable(TestProblem_mysolution ${mysolution_sources})

add_executable(TestProblem_mastersolution ${mastersolution_sources})

add_executable(TestProblem_test_mysolution ${test_sources})

add_executable(TestProblem_test_mastersolution ${test_sources})

target_link_libraries(TestProblem_template
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

target_link_libraries(TestProblem_mysolution
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

target_link_libraries(TestProblem_mastersolution
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

target_link_libraries(TestProblem.solution
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

target_link_libraries(TestProblem.mastersolution
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


target_link_libraries(TestProblem_test_mysolution GTest::main TestProblem.solution)
target_link_libraries(TestProblem_test_mastersolution GTest::main TestProblem.mastersolution)

gtest_discover_tests(TestProblem_test_mysolution)
gtest_discover_tests(TestProblem_test_mastersolution)

