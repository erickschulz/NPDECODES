#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/radauthreetimestepping_main.cc
  ${DIR}/radau_three_timestepping.h
  ${DIR}/radau_three_timestepping.cc
  ${DIR}/radau_three_timestepping_ode.h
  ${DIR}/radau_three_timestepping_ode.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.uscalfe
)
