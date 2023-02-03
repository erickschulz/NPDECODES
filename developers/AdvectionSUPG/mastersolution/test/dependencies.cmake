set(SOURCES
${DIR}/test/advectionsupg_test.cc)
set(LIBRARIES Eigen3::Eigen GTest::gtest_main LF::lf.base LF::lf.mesh LF::lf.mesh.test_utils LF::lf.quad LF::lf.assemble)
