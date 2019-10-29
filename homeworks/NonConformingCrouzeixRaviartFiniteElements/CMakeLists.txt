include(GoogleTest)


set(mastersolution_sources
  mastersolution/compute_cr_l2_error.h
  mastersolution/cr_fe_space.h
  mastersolution/cr_reference_finite_element.cc
  mastersolution/cr_reference_finite_element.h
  mastersolution/cr_types.h
  mastersolution/l2_error_cr_discretization_dirichlet_bvp.h
  mastersolution/non_conforming_crouzeix_raviart_finite_elements_main.cc
  mastersolution/solve_cr_dirichlet_bvp.h
  mastersolution/solve_cr_neumann_bvp.h
)

set(mysolution_sources
  mysolution/compute_cr_l2_error.h
  mysolution/cr_fe_space.h
  mysolution/cr_reference_finite_element.cc
  mysolution/cr_reference_finite_element.h
  mysolution/cr_types.h
  mysolution/l2_error_cr_discretization_dirichlet_bvp.h
  mysolution/non_conforming_crouzeix_raviart_finite_elements_main.cc
  mysolution/solve_cr_dirichlet_bvp.h
  mysolution/solve_cr_neumann_bvp.h
)

set(template_sources
  template/compute_cr_l2_error.h
  template/cr_fe_space.h
  template/cr_reference_finite_element.cc
  template/cr_reference_finite_element.h
  template/cr_types.h
  template/l2_error_cr_discretization_dirichlet_bvp.h
  template/non_conforming_crouzeix_raviart_finite_elements_main.cc
  template/solve_cr_dirichlet_bvp.h
  template/solve_cr_neumann_bvp.h
)

set(mastersolution_test_sources
  test/mastersolution_tests.cc
)

set(mysolution_test_sources
  test/mysolution_tests.cc
)

add_executable(NonConformingCrouzeixRaviart_mastersolution
  ${mastersolution_sources}
)

add_executable(NonConformingCrouzeixRaviart_mysolution
  ${mysolution_sources}
)

add_executable(NonConformingCrouzeixRaviart_template
  ${template_sources}
)

add_library(NonConformingCrouzeixRaviart_mastersolution_lib
  ${mastersolution_sources}
)

add_library(NonConformingCrouzeixRaviart_mysolution_lib
  ${mysolution_sources}
)

add_executable(NonConformingCrouzeixRaviart_test_mastersolution
  ${mastersolution_test_sources}
)

add_executable(NonConformingCrouzeixRaviart_test_mysolution
  ${mysolution_test_sources}
)

target_link_libraries(NonConformingCrouzeixRaviart_mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(NonConformingCrouzeixRaviart_mysolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(NonConformingCrouzeixRaviart_template
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(NonConformingCrouzeixRaviart_mastersolution_lib
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(NonConformingCrouzeixRaviart_mysolution_lib
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

target_link_libraries(NonConformingCrouzeixRaviart_test_mastersolution
  GTest::main
  NonConformingCrouzeixRaviart_mastersolution_lib
)

target_link_libraries(NonConformingCrouzeixRaviart_test_mysolution
  GTest::main
  NonConformingCrouzeixRaviart_mysolution_lib
)


gtest_discover_tests(NonConformingCrouzeixRaviart_test_mastersolution)
gtest_discover_tests(NonConformingCrouzeixRaviart_test_mysolution)
