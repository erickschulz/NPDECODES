# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/lfppdofhandling_main.cc
  ${DIR}/lfppdofhandling.h
  ${DIR}/lfppdofhandling.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
)
