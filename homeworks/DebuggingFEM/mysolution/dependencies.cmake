# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/main.cc
  ${DIR}/local_laplace_qfe.h
  ${DIR}/local_laplace_qfe.cc
  ${DIR}/qfe_interpolator.h
  ${DIR}/qfe_interpolator.cc
  ${DIR}/qfe_provider_tester.h
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)
