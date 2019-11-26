# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/unstablebvp_test.cc
)

set(LIBRARIES
  GTest::gtest_main
  LF::lf.mesh
  LF::lf.refinement
)
