
include(GoogleTest)

set(template_sources
  templates/simple_linear_finite_elements.cc
  templates/simple_linear_finite_elements.h
  templates/tria_mesh_2D.cc
  templates/local_computations.cc
  templates/main.cc
)

set(mysolution_sources
  mysolution/simple_linear_finite_elements.cc
  mysolution/simple_linear_finite_elements.h
  mysolution/tria_mesh_2D.cc
  mysolution/local_computations.cc
  mysolution/main.cc
)

set(mastersolution_sources
  mastersolution/simple_linear_finite_elements.cc
  mastersolution/simple_linear_finite_elements.h
  mastersolution/tria_mesh_2D.cc
  mastersolution/local_computations.cc
  mastersolution/main.cc
)

set(test_sources
  test/simple_linear_finite_elements_tests.cc
)

add_library(SimpleLinearFiniteElements.solution ${mysolution_sources})

add_library(SimpleLinearFiniteElements.mastersolution ${mastersolution_sources})

add_executable(SimpleLinearFiniteElements_template ${template_sources})

add_executable(SimpleLinearFiniteElements_mysolution ${mysolution_sources})

add_executable(SimpleLinearFiniteElements_mastersolution ${mastersolution_sources})

add_executable(SimpleLinearFiniteElements_test_mysolution ${test_sources})
add_executable(SimpleLinearFiniteElements_test_mastersolution ${test_sources})


target_link_libraries(SimpleLinearFiniteElements_template
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
)

target_link_libraries(SimpleLinearFiniteElements_mysolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
)

target_link_libraries(SimpleLinearFiniteElements_mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
)

target_link_libraries(SimpleLinearFiniteElements.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
)

target_link_libraries(SimpleLinearFiniteElements.mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options MathGL
)

target_link_libraries(SimpleLinearFiniteElements_test_mysolution GTest::main SimpleLinearFiniteElements.solution)
target_link_libraries(SimpleLinearFiniteElements_test_mastersolution GTest::main SimpleLinearFiniteElements.mastersolution)


gtest_discover_tests(SimpleLinearFiniteElements_test_mysolution)
gtest_discover_tests(SimpleLinearFiniteElements_test_mastersolution)




