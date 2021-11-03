# Dependencies of mastersolution tests:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/quasiinterpolation_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.base
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.uscalfe
)
