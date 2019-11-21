#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/symplectictimesteppingwaves_main.cc
  ${DIR}/symplectictimesteppingwaves.h
  ${DIR}/symplectictimesteppingwaves.cc
  ${DIR}/symplectictimesteppingwaves_assemble.h
  ${DIR}/symplectictimesteppingwaves_ode.h
  ${DIR}/symplectictimesteppingwaves_ode.cc
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
