#if SOLUTION
# Dependencies of mastersolution tests:
#else
# Dependencies of mysolution tests:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/discontinuousgalerkin1d_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
)
