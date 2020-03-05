# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/cr_l2_error.h
  ${DIR}/cr_fe_space.h
  ${DIR}/nonconformingcrouzeixraviartfiniteelements.cc
  ${DIR}/nonconformingcrouzeixraviartfiniteelements.h
  ${DIR}/l2_error_cr_discretization_dirichlet_bvp.h
  ${DIR}/nonconformingcrouzeixraviartfiniteelements_main.cc
  ${DIR}/solve_cr_dirichlet_bvp.h
  ${DIR}/solve_cr_neumann_bvp.h
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.io
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.uscalfe
)
