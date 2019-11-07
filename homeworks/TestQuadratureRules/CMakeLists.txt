include(GoogleTest)

set(mastersolution_sources
  mastersolution/test_quad_rules.cc
  mastersolution/test_quad_rules.h
  mastersolution/TestQuadratureRules.cc
)


set(mastersolution_test_sources
  test/test_quad_rules_test_mastersolution.cc
)


# Libraries
add_library(TestQuadratureRules_mastersolution.solution ${mastersolution_sources})



# Executables
add_executable(TestQuadratureRules_mastersolution ${mastersolution_sources})

add_executable(TestQuadratureRules_test_mastersolution ${mastersolution_test_sources})


target_link_libraries(TestQuadratureRules_mastersolution
	PUBLIC Eigen3::Eigen
	LF::lf.base
	LF::lf.quad
)


target_link_libraries(TestQuadratureRules_mastersolution.solution
	PUBLIC Eigen3::Eigen
	LF::lf.base
	LF::lf.quad
)


target_link_libraries(TestQuadratureRules_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::gtest_main
	TestQuadratureRules_mastersolution.solution
	LF::lf.base
	LF::lf.quad
)

gtest_discover_tests(TestQuadratureRules_test_mastersolution)
