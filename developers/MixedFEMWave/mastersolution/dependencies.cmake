#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/mixedfemwave_main.cc
  ${DIR}/mixedfemwave.h
  ${DIR}/mixedfemwave.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.base
  LF::lf.io
  LF::lf.geometry
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.uscalfe
)
