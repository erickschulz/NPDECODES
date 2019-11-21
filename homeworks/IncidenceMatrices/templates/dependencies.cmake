# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/main.cc
  ${DIR}/incidence_mat.h
  ${DIR}/incidence_mat.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.mesh.hybrid2d
)
