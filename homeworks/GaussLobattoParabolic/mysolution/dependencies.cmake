# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/gausslobattoparabolic_main.cc
  ${DIR}/gausslobattoparabolic.h
  ${DIR}/gausslobattoparabolic.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.uscalfe
)
