# Dependencies of mastersolution tests:

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/rk3prey_test.cc
  ${DIR}/rk3prey.h
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
)
