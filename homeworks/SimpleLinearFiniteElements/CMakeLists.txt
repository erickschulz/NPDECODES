
include(GoogleTest)


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


add_library(SimpleLinearFiniteElements.mastersolution ${mastersolution_sources})


add_executable(SimpleLinearFiniteElements_mastersolution ${mastersolution_sources})

add_executable(SimpleLinearFiniteElements_test_mastersolution ${test_sources})


target_link_libraries(SimpleLinearFiniteElements_mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options 
)


target_link_libraries(SimpleLinearFiniteElements.mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options 
)

target_link_libraries(SimpleLinearFiniteElements_test_mastersolution GTest::gtest_main SimpleLinearFiniteElements.mastersolution)


gtest_discover_tests(SimpleLinearFiniteElements_test_mastersolution)




