include(GoogleTest)

set(template_sources
  templates/maximum_principle_main.cc
  templates/maximum_principle.h
  templates/maximum_principle.cc
)

set(mastersolution_sources
  mastersolution/maximum_principle_main.cc
  mastersolution/maximum_principle.h
  mastersolution/maximum_principle.cc
)

set(mysolution_sources
  mysolution/maximum_principle_main.cc
  mysolution/maximum_principle.h
  mysolution/maximum_principle.cc
)

set(mastersolution_test_sources
  test/maximum_principle_test_mastersolution.cc
)

set(mysolution_test_sources
  test/maximum_principle_test_mysolution.cc
)

# Libraries
add_library(MaximumPrinciple_mastersolution.solution ${mastersolution_sources})

add_library(MaximumPrinciple_mysolution.solution ${mysolution_sources})

# Executables
add_executable(MaximumPrinciple_template ${template_sources})

add_executable(MaximumPrinciple_mysolution ${mysolution_sources})

add_executable(MaximumPrinciple_mastersolution ${mastersolution_sources})

add_executable(MaximumPrinciple_test_mastersolution ${mastersolution_test_sources})

add_executable(MaximumPrinciple_test_mysolution ${mysolution_test_sources})


target_link_libraries(MaximumPrinciple_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(MaximumPrinciple_mysolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(MaximumPrinciple_template
	PUBLIC Eigen3::Eigen
)

target_link_libraries(MaximumPrinciple_mastersolution.solution
  PUBLIC Eigen3::Eigen 
  Boost::system
)

target_link_libraries(MaximumPrinciple_mysolution.solution
  PUBLIC Eigen3::Eigen 
  Boost::system
)

target_link_libraries(MaximumPrinciple_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	MaximumPrinciple_mastersolution.solution
)

target_link_libraries(MaximumPrinciple_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	MaximumPrinciple_mysolution.solution
)

gtest_discover_tests(MaximumPrinciple_test_mysolution)
gtest_discover_tests(MaximumPrinciple_test_mastersolution)
