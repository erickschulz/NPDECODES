# Dependencies of mastersolution tests:

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/projectionontogradients_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.base
  LF::lf.mesh.hybrid2d
  LF::lf.mesh
  LF::lf.mesh.test_utils
  LF::lf.uscalfe
)
