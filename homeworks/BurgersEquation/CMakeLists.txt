include(GoogleTest)

set(mastersolution_sources
  mastersolution/burgersequation_main.cc
  mastersolution/burgersequation.h
  mastersolution/burgersequation.cc
)
configure_file(mastersolution/plot_solution.py plot_solution_mastersolution.py COPYONLY)
configure_file(mastersolution/plot_error.py plot_error_mastersolution.py COPYONLY)

set(mastersolution_test_sources
  test/burgersequation_test_mastersolution.cc
)

# Executables
add_executable(BurgersEquation_mastersolution ${mastersolution_sources})
add_executable(BurgersEquation_test_mastersolution ${mastersolution_test_sources})

# Libraries
add_library(BurgersEquation_mastersolution.solution ${mastersolution_sources})

# Linking
target_link_libraries(BurgersEquation_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(BurgersEquation_mastersolution.solution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(BurgersEquation_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::gtest_main
	BurgersEquation_mastersolution.solution
)

gtest_discover_tests(BurgersEquation_test_mastersolution)
