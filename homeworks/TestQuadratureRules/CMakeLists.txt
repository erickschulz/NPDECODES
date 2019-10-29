include(GoogleTest)

set(mastersolution_sources
  mastersolution/test_quad_rules.cc
  mastersolution/test_quad_rules.h
  mastersolution/TestQuadratureRules.cc
)

set(mysolution_sources
  mysolution/test_quad_rules.cc
  mysolution/test_quad_rules.h
  mysolution/TestQuadratureRules.cc
)

set(template_sources
  templates/test_quad_rules.cc
  templates/test_quad_rules.h
  templates/TestQuadratureRules.cc
)

set(mastersolution_test_sources
  test/test_quad_rules_test_mastersolution.cc
)

set(mysolution_test_sources
  test/test_quad_rules_test_mysolution.cc
)

# Libraries
add_library(TestQuadratureRules_mastersolution.solution ${mastersolution_sources})

add_library(TestQuadratureRules_mysolution.solution ${mysolution_sources})


# Executables
add_executable(TestQuadratureRules_template ${template_sources})

add_executable(TestQuadratureRules_mysolution ${mysolution_sources})

add_executable(TestQuadratureRules_mastersolution ${mastersolution_sources})

add_executable(TestQuadratureRules_test_mastersolution ${mastersolution_test_sources})

add_executable(TestQuadratureRules_test_mysolution ${mysolution_test_sources})

target_link_libraries(TestQuadratureRules_mastersolution
	PUBLIC Eigen3::Eigen
	Boost::system
	LF::lf.base
	LF::lf.quad
)

target_link_libraries(TestQuadratureRules_mysolution
	PUBLIC Eigen3::Eigen
	Boost::system
	LF::lf.base
	LF::lf.quad
)

target_link_libraries(TestQuadratureRules_template
	PUBLIC Eigen3::Eigen
	Boost::system
	LF::lf.base
	LF::lf.quad
)

target_link_libraries(TestQuadratureRules_mastersolution.solution
	PUBLIC Eigen3::Eigen
	Boost::system
	LF::lf.base
	LF::lf.quad
)

target_link_libraries(TestQuadratureRules_mysolution.solution
	PUBLIC Eigen3::Eigen
	Boost::system
	LF::lf.base
	LF::lf.quad
)

target_link_libraries(TestQuadratureRules_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	TestQuadratureRules_mastersolution.solution
	Boost::system
	LF::lf.base
	LF::lf.quad
)

target_link_libraries(TestQuadratureRules_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	TestQuadratureRules_mysolution.solution
	Boost::system
	LF::lf.base
	LF::lf.quad
)

gtest_discover_tests(TestQuadratureRules_test_mysolution)
gtest_discover_tests(TestQuadratureRules_test_mastersolution)
