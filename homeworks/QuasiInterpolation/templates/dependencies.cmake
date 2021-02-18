# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/quasiinterpolation_main.cc
  ${DIR}/quasiinterpolation.h
  ${DIR}/quasiinterpolation.cc
  ${DIR}/iohelper.h
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.base
  LF::lf.io
  LF::lf.quad
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.uscalfe
)
