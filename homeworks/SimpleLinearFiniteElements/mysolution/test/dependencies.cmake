# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/simple_linear_finite_elements_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
)
