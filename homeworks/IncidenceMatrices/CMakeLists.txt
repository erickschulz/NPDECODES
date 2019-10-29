include(GoogleTest)

set(template_sources
  templates/main.cc
  templates/incidence_mat.cc
  templates/incidence_mat.h
)

set(mysolution_sources
  mysolution/main.cc
  mysolution/incidence_mat.cc
  mysolution/incidence_mat.h
)

set(mastersolution_sources
  mastersolution/main.cc
  mastersolution/incidence_mat.cc
  mastersolution/incidence_mat.h
)

set(test_mysolution_sources
  test/test_incidence_mat.cc
  mysolution/incidence_mat.cc
  mysolution/incidence_mat.h
)

set(test_mastersolution_sources
  test/test_incidence_mat.cc
  mastersolution/incidence_mat.cc
  mastersolution/incidence_mat.h
)

add_executable(IncidenceMatrices_template ${template_sources})
add_executable(IncidenceMatrices_mysolution ${mysolution_sources})
add_executable(IncidenceMatrices_mastersolution ${mastersolution_sources})
add_executable(IncidenceMatrices_test_mysolution ${test_mysolution_sources})
add_executable(IncidenceMatrices_test_mastersolution ${test_mastersolution_sources})

target_link_libraries(IncidenceMatrices_template
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
)

target_link_libraries(IncidenceMatrices_mysolution
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
)

target_link_libraries(IncidenceMatrices_mastersolution
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
)

target_link_libraries(IncidenceMatrices_test_mysolution
  GTest::main
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
)

target_link_libraries(IncidenceMatrices_test_mastersolution
  GTest::main
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.hybrid2d
)

gtest_discover_tests(IncidenceMatrices_test_mysolution)
gtest_discover_tests(IncidenceMatrices_test_mastersolution)

add_custom_target(IncidenceMatrices_maketemplate
    COMMAND ${CMAKE_SOURCE_DIR}/strip-comments.sh
)
