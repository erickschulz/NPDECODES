# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/advectionfv2d_main.cc
  ${DIR}/advectionfv2d.h
  ${DIR}/advectionfv2d.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.refinement
)
