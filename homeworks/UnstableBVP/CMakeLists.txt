include(GoogleTest)

set(libraries
  Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.uscalfe
)

set(mastersolution_sources
  mastersolution/unstablebvp.h
  mastersolution/unstablebvp.cc
  mastersolution/unstablebvp_main.cc
)

set(test_mastersolution_sources
  mastersolution/unstablebvp.h
  mastersolution/unstablebvp.cc
  test/unstablebvp_test_mastersolution.cc
)

add_executable(UnstableBVP_mastersolution ${mastersolution_sources})

target_link_libraries(UnstableBVP_mastersolution PUBLIC ${libraries})
target_compile_features(UnstableBVP_mastersolution PUBLIC cxx_std_17)

add_executable(test_mastersolution_UnstableBVP ${test_mastersolution_sources})

target_link_libraries(test_mastersolution_UnstableBVP PUBLIC ${libraries} GTest::gtest_main)
target_compile_features(test_mastersolution_UnstableBVP PUBLIC cxx_std_17)

gtest_discover_tests(test_mastersolution_UnstableBVP)
