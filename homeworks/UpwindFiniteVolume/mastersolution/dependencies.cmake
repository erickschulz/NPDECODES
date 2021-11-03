# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/upwindfinitevolume_main.cc
  ${DIR}/upwindfinitevolume.h
  ${DIR}/upwindfinitevolume.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
)
