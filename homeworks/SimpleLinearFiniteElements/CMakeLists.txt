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

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

add_executable(SimpleLinearFiniteElements_mastersolution ${mastersolution_sources})
add_executable(SimpleLinearFiniteElements_test_mastersolution ${test_sources})

add_library(SimpleLinearFiniteElements.mastersolution ${mastersolution_sources})

target_link_libraries(SimpleLinearFiniteElements_mastersolution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(SimpleLinearFiniteElements.mastersolution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(SimpleLinearFiniteElements_test_mastersolution
  GTest::gtest_main
  SimpleLinearFiniteElements.mastersolution
)

gtest_discover_tests(SimpleLinearFiniteElements_test_mastersolution)
