
include(GoogleTest)

set(mastersolution_sources
  mastersolution/ansiotropic_diffusion_element_matrix_provider.cc
  mastersolution/ansiotropic_diffusion_element_matrix_provider.h
  mastersolution/fe_source_elem_vec_provider.cc
  mastersolution/fe_source_elem_vec_provider.h
  mastersolution/impedance_boundary_edge_matrix_provider.h
  mastersolution/impedance_boundary_edge_matrix_provider.cc
  mastersolution/main.cc
)

set(mastersolution_test
  mastersolution/ansiotropic_diffusion_element_matrix_provider.cc
  mastersolution/ansiotropic_diffusion_element_matrix_provider.h
  mastersolution/fe_source_elem_vec_provider.cc
  mastersolution/fe_source_elem_vec_provider.h
  mastersolution/impedance_boundary_edge_matrix_provider.h
  mastersolution/impedance_boundary_edge_matrix_provider.cc
  mastersolution/parametric_element_matrices_test.cc
)

# set(mysolution_sources
#   mysolution/ansiotropic_diffusion_element_matrix_provider.cc
#   mysolution/ansiotropic_diffusion_element_matrix_provider.h
#   mysolution/fe_source_elem_vec_provider.cc
#   mysolution/fe_source_elem_vec_provider.h
#   mysolution/impedance_boundary_edge_matrix_provider.h
#   mysolution/impedance_boundary_edge_matrix_provider.cc
#   mysolution/main.cc
# )

# set(mysolution_test
#   mysolution/ansiotropic_diffusion_element_matrix_provider.cc
#   mysolution/ansiotropic_diffusion_element_matrix_provider.h
#   mysolution/fe_source_elem_vec_provider.cc
#   mysolution/fe_source_elem_vec_provider.h
#   mysolution/impedance_boundary_edge_matrix_provider.h
#   mysolution/impedance_boundary_edge_matrix_provider.cc
#   mysolution/parametric_element_matrices_test.cc
# )

# set(templates_sources
#   templates/ansiotropic_diffusion_element_matrix_provider.cc
#   templates/ansiotropic_diffusion_element_matrix_provider.h
#   templates/fe_source_elem_vec_provider.cc
#   templates/fe_source_elem_vec_provider.h
#   templates/impedance_boundary_edge_matrix_provider.h
#   templates/impedance_boundary_edge_matrix_provider.cc
#   templates/main.cc
# )

# set(templates_test
#   templates/ansiotropic_diffusion_element_matrix_provider.cc
#   templates/ansiotropic_diffusion_element_matrix_provider.h
#   templates/fe_source_elem_vec_provider.cc
#   templates/fe_source_elem_vec_provider.h
#   templates/impedance_boundary_edge_matrix_provider.h
#   templates/impedance_boundary_edge_matrix_provider.cc
#   templates/parametric_element_matrices_test.cc
# )

# Libraries
add_executable(ParametricElementMatrices_mastersolution ${mastersolution_sources})

add_executable(ParametricElementMatrices_test_mastersolution ${mastersolution_test})

# add_executable(ParametricElementMatrices_mysolution ${mysolution_sources})

# add_executable(ParametricElementMatrices_test_mysolution ${mysolution_test})

# add_executable(ParametricElementMatrices_templates ${templates_sources})

# add_executable(ParametricElementMatrices_test_templates ${templates_test})

target_link_libraries(ParametricElementMatrices_mastersolution
	PUBLIC Eigen3::Eigen
	Boost::system
	LF::lf.mesh
	LF::lf.mesh.test_utils
	LF::lf.assemble
	LF::lf.uscalfe
	LF::lf.quad
	LF::lf.geometry
)

target_link_libraries(ParametricElementMatrices_test_mastersolution
	PUBLIC Eigen3::Eigen 
	Boost::system
	LF::lf.mesh
	LF::lf.mesh.test_utils
	LF::lf.assemble
	LF::lf.uscalfe
	LF::lf.quad
	LF::lf.geometry
)

# target_link_libraries(ParametricElementMatrices_mysolution
# 	PUBLIC Eigen3::Eigen
# 	Boost::system
# 	LF::lf.mesh
# 	LF::lf.mesh.test_utils
# 	LF::lf.assemble
# 	LF::lf.uscalfe
# 	LF::lf.quad
# 	LF::lf.geometry
# )

# target_link_libraries(ParametricElementMatrices_test_mysolution
# 	PUBLIC Eigen3::Eigen 
# 	Boost::system
# 	LF::lf.mesh
# 	LF::lf.mesh.test_utils
# 	LF::lf.assemble
# 	LF::lf.uscalfe
# 	LF::lf.quad
# 	LF::lf.geometry
# )

# target_link_libraries(ParametricElementMatrices_templates
# 	PUBLIC Eigen3::Eigen
# 	Boost::system
# 	LF::lf.mesh
# 	LF::lf.mesh.test_utils
# 	LF::lf.assemble
# 	LF::lf.uscalfe
# 	LF::lf.quad
# 	LF::lf.geometry
# )

# target_link_libraries(ParametricElementMatrices_test_templates
# 	PUBLIC Eigen3::Eigen 
# 	Boost::system
# 	LF::lf.mesh
# 	LF::lf.mesh.test_utils
# 	LF::lf.assemble
# 	LF::lf.uscalfe
# 	LF::lf.quad
# 	LF::lf.geometry
# )


gtest_discover_tests(ParametricElementMatrices_test_mastersolution)
# gtest_discover_tests(ParametricElementMatrices_test_mysolution)
# gtest_discover_tests(ParametricElementMatrices_test_templates)


