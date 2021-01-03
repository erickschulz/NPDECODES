set(SOURCES
  ${DIR}/test/expfittedupwind_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.uscalfe
)
