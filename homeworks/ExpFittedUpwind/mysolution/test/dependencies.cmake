set(SOURCES
  ${DIR}/test/expfittedupwind_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.base
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.uscalfe
)
