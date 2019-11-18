# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/boundarylength_main.cc
  ${DIR}/boundarylength.h
  ${DIR}/boundarylength.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.utils
  LF::lf.mesh.hybrid2d
)
