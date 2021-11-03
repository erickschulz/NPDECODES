# Dependencies of mastersolution tests:

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/gradientflow_test.cc
  ${DIR}/gradientflow.h
  ${DIR}/gradientflow.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
)
