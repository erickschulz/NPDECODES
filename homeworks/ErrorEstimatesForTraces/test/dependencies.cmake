# Dependencies of mastersolution tests:

# DIR will be provided by the calling file.

set(SOURCES
  test/tee_lapl_robin_assembly_test_${DIR}.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.mesh.utils
  LF::lf.mesh.hybrid2d
  LF::lf.io
)
