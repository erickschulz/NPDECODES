# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/sdirkmethodoflines_main.cc
  ${DIR}/sdirkmethodoflines.h
  ${DIR}/sdirkmethodoflines.cc
  ${DIR}/sdirkmethodoflines_ode.h
  ${DIR}/sdirkmethodoflines_ode.cc
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
