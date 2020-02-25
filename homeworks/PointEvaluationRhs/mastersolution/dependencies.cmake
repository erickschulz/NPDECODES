# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/pointevaluationrhs_main.cc
  ${DIR}/pointevaluationrhs.h
  ${DIR}/pointevaluationrhs.cc
  ${DIR}/pointevaluationrhs_norms.h
  ${DIR}/pointevaluationrhs_norms.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.uscalfe
)
