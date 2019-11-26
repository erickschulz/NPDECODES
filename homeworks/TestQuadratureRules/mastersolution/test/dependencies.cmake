# Dependencies of mastersolution tests:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/test_quad_rules_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.base
  LF::lf.quad
)
