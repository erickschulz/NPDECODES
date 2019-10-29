include(GoogleTest)

set(template_sources
templates/PointEvaluationRhs_main.cc
templates/norms.cc
templates/pointEvaluation.cc
)
set(mysolution_sources
mysolution/PointEvaluationRhs_main.cc
mysolution/norms.cc
mysolution/pointEvaluation.cc
)
set(mastersolution_sources
mastersolution/norms.cc
mastersolution/pointEvaluation.cc
mastersolution/PointEvaluationRhs_main.cc
)
set(master_test_sources
        mastersolution/PointEvaluationRhs_test.cc
)
set(my_test_sources
  mysolution/PointEvaluationRhs_test.cc
)
set(template_test_sources
  templates/PointEvaluationRhs_test.cc
)

add_executable(PointEvaluationRhs_template ${template_sources})
add_executable(PointEvaluationRhs_mysolution ${mysolution_sources})
add_executable(PointEvaluationRhs_mastersolution ${mastersolution_sources})

target_link_libraries(PointEvaluationRhs_template
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)
target_link_libraries(PointEvaluationRhs_mysolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(PointEvaluationRhs_mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
  LF::lf.mesh
  LF::lf.mesh.utils
        LF::lf.refinement
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

#Unit tests
add_executable(PointEvaluationRhs_test_mysolution ${my_test_sources})
add_executable(PointEvaluationRhs_test_mastersolution ${master_test_sources})
 add_executable(PointEvaluationRhs_test_template ${template_test_sources})

add_library(PointEvaluationRhs_my.solution ${mysolution_sources})
add_library(PointEvaluationRhs_master.solution ${mastersolution_sources})
add_library(PointEvaluationRhs_template.solution ${template_sources})

target_link_libraries(PointEvaluationRhs_my.solution
 PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
 LF::lf.mesh
 LF::lf.mesh.utils
 LF::lf.mesh.test_utils
       LF::lf.refinement
  LF::lf.mesh.hybrid2d
 LF::lf.assemble
 LF::lf.io
 LF::lf.uscalfe
 )
target_link_libraries(PointEvaluationRhs_master.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
        LF::lf.refinement
  LF::lf.io
  LF::lf.uscalfe
  )
target_link_libraries(PointEvaluationRhs_template.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  )

target_link_libraries(PointEvaluationRhs_test_mysolution
GTest::main
PointEvaluationRhs_my.solution
)
target_link_libraries(PointEvaluationRhs_test_mastersolution
GTest::main 
PointEvaluationRhs_master.solution
)
target_link_libraries(PointEvaluationRhs_test_template
GTest::main
PointEvaluationRhs_template.solution
)

gtest_discover_tests(PointEvaluationRhs_test_mastersolution)
gtest_discover_tests(PointEvaluationRhs_test_mysolution)
gtest_discover_tests(PointEvaluationRhs_test_template)


