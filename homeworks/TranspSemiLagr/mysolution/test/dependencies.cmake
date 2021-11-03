set(SOURCES
  ${DIR}/test/transpsemilagr_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.mesh
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.uscalfe
)