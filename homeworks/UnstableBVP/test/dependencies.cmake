# Dependencies of mastersolution tests:

# DIR will be provided by the calling file.

set(SOURCES
  test/unstablebvp_test_${DIR}.cc
)

set(LIBRARIES
  GTest::gtest_main
  LF::lf.mesh
  LF::lf.refinement
)
