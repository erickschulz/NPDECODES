# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/lfppdofhandling_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
)
