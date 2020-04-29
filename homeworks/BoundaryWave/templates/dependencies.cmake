# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/boundarywave_main.cc
  ${DIR}/boundarywave.h
  ${DIR}/boundarywave.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.uscalfe
)
