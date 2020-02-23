# Dependencies of mastersolution tests:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/PointEvaluationRhs_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
)
