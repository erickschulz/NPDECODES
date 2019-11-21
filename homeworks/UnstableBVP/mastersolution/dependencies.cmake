# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/unstablebvp_main.cc
  ${DIR}/unstablebvp.h
  ${DIR}/unstablebvp.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.uscalfe
)
