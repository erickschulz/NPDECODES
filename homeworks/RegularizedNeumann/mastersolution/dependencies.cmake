# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/regularizedneumann_main.cc
  ${DIR}/regularizedneumann_get_galerkin_lse.h
  ${DIR}/regularizedneumann.h
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.mesh
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.uscalfe
)
