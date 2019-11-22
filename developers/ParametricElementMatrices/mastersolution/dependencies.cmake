#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/parametricelementmatrices_main.cc
  ${DIR}/ansiotropic_diffusion_element_matrix_provider.h
  ${DIR}/ansiotropic_diffusion_element_matrix_provider.cc
  ${DIR}/fe_source_elem_vec_provider.h
  ${DIR}/fe_source_elem_vec_provider.cc
  ${DIR}/impedance_boundary_edge_matrix_provider.h
  ${DIR}/impedance_boundary_edge_matrix_provider.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh
  LF::lf.mesh.test_utils
  LF::lf.quad
  LF::lf.uscalfe
)
