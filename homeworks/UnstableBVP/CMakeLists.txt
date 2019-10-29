include(GoogleTest)

set(libraries
  Eigen3::Eigen
  Boost::filesystem
  Boost::system
  Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

set(sources
  unstable_bvp.h
  unstable_bvp.cc
  unstable_bvp_main.cc
)

set(test_sources
  unstable_bvp.h
  unstable_bvp.cc
  test/test_unstable_bvp.cc
)

add_executable(UnstableBVP ${sources})

target_link_libraries(UnstableBVP PUBLIC ${libraries})
target_compile_features(UnstableBVP PUBLIC cxx_std_17)

add_executable(test_UnstableBVP ${test_sources})

target_link_libraries(test_UnstableBVP PUBLIC ${libraries} GTest::main)
target_compile_features(test_UnstableBVP PUBLIC cxx_std_17)

gtest_discover_tests(test_UnstableBVP)
