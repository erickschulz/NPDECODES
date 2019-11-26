#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/electrostaticforce_main.cc
  ${DIR}/electrostaticforce.h
  ${DIR}/electrostaticforce.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.uscalfe
)
