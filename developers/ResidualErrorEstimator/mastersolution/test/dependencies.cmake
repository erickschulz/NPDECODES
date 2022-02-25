set(SOURCES
  ${DIR}/test/residualerrorestimator_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.base
  LF::lf.geometry
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.quad
  LF::lf.io
  LF::lf.fe
  LF::lf.uscalfe
)

