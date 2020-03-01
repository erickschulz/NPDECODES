#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/handlingdofs_main.cc
  ${DIR}/handlingdofs.h
  ${DIR}/handlingdofs.cc
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
