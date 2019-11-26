#if SOLUTION
# Dependencies of mastersolution tests:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/boundarylength_test.cc
)

set(LIBRARIES
  GTest::gtest_main
  LF::lf.mesh.test_utils
)
