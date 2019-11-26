#if SOLUTION
# Dependencies of mastersolution tests:
#else
# Add your custom dependencies here:
#endif

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  test/clempiricflux_test_${DIR}.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
)
