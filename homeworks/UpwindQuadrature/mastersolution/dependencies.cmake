# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/upwindquadrature_main.cc
  ${DIR}/upwindquadrature.h
  ${DIR}/upwindquadrature.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.uscalfe
)
