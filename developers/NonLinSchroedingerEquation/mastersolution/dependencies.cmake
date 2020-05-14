#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/nonlinschroedingerequation_main.cc
  ${DIR}/nonlinschroedingerequation.h
  ${DIR}/nonlinschroedingerequation.cc
  ${DIR}/propagator.h
  ${DIR}/propagator.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.uscalfe
)
