# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/convblfmatrixprovider.h
  ${DIR}/convblfmatrixprovider.cc
  ${DIR}/convblfmatrixprovider_main.cc
  )

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.io
  LF::lf.geometry
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.fe
  LF::lf.uscalfe
)
