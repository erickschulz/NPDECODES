# Dependencies of mastersolution tests:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/linfereactdiff_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.base
  LF::lf.io
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
)
