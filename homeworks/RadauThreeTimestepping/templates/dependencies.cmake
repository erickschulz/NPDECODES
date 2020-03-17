# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/radauthreetimestepping_main.cc
  ${DIR}/radauthreetimestepping.h
  ${DIR}/radauthreetimestepping.cc
  ${DIR}/radauthreetimesteppingode.h
  ${DIR}/radauthreetimesteppingode.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.uscalfe
)
