include(GoogleTest)



set(mastersolution_sources
  mastersolution/main.cc
  mastersolution/incidence_mat.cc
  mastersolution/incidence_mat.h
)


set(test_mastersolution_sources
  test/test_master_incidence_mat.cc
  mastersolution/incidence_mat.cc
  mastersolution/incidence_mat.h
)


add_executable(IncidenceMatrices_mastersolution ${mastersolution_sources})
add_executable(IncidenceMatrices_test_mastersolution ${test_mastersolution_sources})


target_link_libraries(IncidenceMatrices_mastersolution
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
)



target_link_libraries(IncidenceMatrices_test_mastersolution
  GTest::gtest_main
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
)

gtest_discover_tests(IncidenceMatrices_test_mastersolution)


