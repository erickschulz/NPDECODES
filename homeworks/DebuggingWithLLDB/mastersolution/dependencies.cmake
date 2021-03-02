# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/debuggingwithlldb_main.cc
  ${DIR}/debuggingwithlldb.cc
  ${DIR}/debuggingwithlldb.h
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.io
 )
