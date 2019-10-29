include(GoogleTest)

set(template_sources
  templates/burgersequation_main.cc
  templates/burgersequation.h
  templates/burgersequation.cc
)

set(mastersolution_sources
  mastersolution/burgersequation_main.cc
  mastersolution/burgersequation.h
  mastersolution/burgersequation.cc
)
configure_file(mastersolution/plot_solution.py plot_solution_mastersolution.py COPYONLY)
configure_file(mastersolution/plot_error.py plot_error_mastersolution.py COPYONLY)

set(mysolution_sources
  mysolution/burgersequation_main.cc
  mysolution/burgersequation.h
  mysolution/burgersequation.cc
)

set(mastersolution_test_sources
  test/burgersequation_test_mastersolution.cc
)

set(mysolution_test_sources
  test/burgersequation_test_mysolution.cc
)

# Libraries
add_library(BurgersEquation_mastersolution.solution ${mastersolution_sources})

add_library(BurgersEquation_mysolution.solution ${mysolution_sources})

# Executables
add_executable(BurgersEquation_template ${template_sources})

add_executable(BurgersEquation_mysolution ${mysolution_sources})

add_executable(BurgersEquation_mastersolution ${mastersolution_sources})

add_executable(BurgersEquation_test_mastersolution ${mastersolution_test_sources})

add_executable(BurgersEquation_test_mysolution ${mysolution_test_sources})


target_link_libraries(BurgersEquation_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(BurgersEquation_mysolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(BurgersEquation_template
	PUBLIC Eigen3::Eigen
)

target_link_libraries(BurgersEquation_mastersolution.solution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(BurgersEquation_mysolution.solution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(BurgersEquation_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	BurgersEquation_mastersolution.solution
)

target_link_libraries(BurgersEquation_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	BurgersEquation_mysolution.solution
)

gtest_discover_tests(BurgersEquation_test_mysolution)
gtest_discover_tests(BurgersEquation_test_mastersolution)
