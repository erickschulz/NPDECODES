set(SOURCES
${DIR}/advectionsupg_main.cc
${DIR}/advectionsupg.h)
set(LIBRARIES Eigen3::Eigen LF::lf.base LF::lf.mesh LF::lf.mesh.test_utils LF::lf.quad LF::lf.assemble)
