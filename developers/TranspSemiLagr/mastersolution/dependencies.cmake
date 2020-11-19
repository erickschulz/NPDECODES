#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

set(SOURCES
  ${DIR}/local_assembly.h
  ${DIR}/transpsemilagr.h
  ${DIR}/transpsemilagr_main.cc
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
  LF::lf.uscalfe
)
