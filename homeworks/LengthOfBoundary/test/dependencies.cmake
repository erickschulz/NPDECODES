#if SOLUTION
# Dependencies of mastersolution tests:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  test/boundarylength_test_${DIR}.cc
)

set(LIBRARIES
  GTest::gtest_main
  LF::lf.mesh.test_utils
)
