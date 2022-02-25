# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/stableevaluationatapoint_main.cc
  ${DIR}/stableevaluationatapoint.h
  ${DIR}/stableevaluationatapoint.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.fe
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.quad
  LF::lf.uscalfe
)
