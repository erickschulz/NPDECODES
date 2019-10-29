include(GoogleTest)

set(template_sources
templates/gradprojection.h
templates/main.cc
)
set(mysolution_sources
mysolution/gradprojection.h
mysolution/main.cc
)
set(mastersolution_sources
mastersolution/gradprojection.h
mastersolution/main.cc
)
set(master_test_sources
  mastersolution/test.cc
)
set(my_test_sources
  mysolution/test.cc
)
set(template_test_sources
  templates/test.cc
)

add_executable(ProjectionOntoGradients_template ${template_sources})
add_executable(ProjectionOntoGradients_mysolution ${mysolution_sources})
add_executable(ProjectionOntoGradients_mastersolution ${mastersolution_sources})

target_link_libraries(ProjectionOntoGradients_template
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)
target_link_libraries(ProjectionOntoGradients_mysolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.io
  LF::lf.uscalfe
)
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
add_executable(ProjectionOntoGradients_test_mysolution ${my_test_sources})
add_executable(ProjectionOntoGradients_test_mastersolution ${master_test_sources})
add_executable(ProjectionOntoGradients_test_template ${template_test_sources})

add_library(ProjectionOntoGradients_my.solution ${mysolution_sources})
add_library(ProjectionOntoGradients_master.solution ${mastersolution_sources})
add_library(ProjectionOntoGradients_template.solution ${template_sources})

target_link_libraries(ProjectionOntoGradients_my.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  )
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
target_link_libraries(ProjectionOntoGradients_template.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  )

target_link_libraries(ProjectionOntoGradients_test_mysolution
GTest::main 
ProjectionOntoGradients_my.solution
)
target_link_libraries(ProjectionOntoGradients_test_mastersolution
GTest::main 
ProjectionOntoGradients_master.solution
)
target_link_libraries(ProjectionOntoGradients_test_template
GTest::main 
ProjectionOntoGradients_template.solution
)

gtest_discover_tests(ProjectionOntoGradients_test_mastersolution)
gtest_discover_tests(ProjectionOntoGradients_test_mysolution)
gtest_discover_tests(ProjectionOntoGradients_test_template)


