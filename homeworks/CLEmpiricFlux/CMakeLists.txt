include(GoogleTest)

set(template_sources
  templates/clempiricflux_main.cc
  templates/clempiricflux.h
  templates/clempiricflux.cc
  templates/uniformcubicspline.h
  templates/uniformcubicspline.cc
  templates/solvecauchyproblem.h
  templates/solvecauchyproblem.cc
)

set(mastersolution_sources
  mastersolution/clempiricflux_main.cc
  mastersolution/clempiricflux.h
  mastersolution/clempiricflux.cc
  mastersolution/uniformcubicspline.h
  mastersolution/uniformcubicspline.cc
  mastersolution/solvecauchyproblem.h
  mastersolution/solvecauchyproblem.cc
)

set(mysolution_sources
  mysolution/clempiricflux_main.cc
  mysolution/clempiricflux.h
  mysolution/clempiricflux.cc
  mysolution/uniformcubicspline.h
  mysolution/uniformcubicspline.cc
  mysolution/solvecauchyproblem.h
  mysolution/solvecauchyproblem.cc
)

set(mastersolution_test_sources
  test/clempiricflux_test_mastersolution.cc
)

set(mysolution_test_sources
  test/clempiricflux_test_mysolution.cc
)

# Libraries
add_library(CLEmpiricFlux_mastersolution.solution ${mastersolution_sources})

add_library(CLEmpiricFlux_mysolution.solution ${mysolution_sources})

# Executables
add_executable(CLEmpiricFlux_template ${template_sources})

add_executable(CLEmpiricFlux_mysolution ${mysolution_sources})

add_executable(CLEmpiricFlux_mastersolution ${mastersolution_sources})

add_executable(CLEmpiricFlux_test_mastersolution ${mastersolution_test_sources})

add_executable(CLEmpiricFlux_test_mysolution ${mysolution_test_sources})


target_link_libraries(CLEmpiricFlux_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(CLEmpiricFlux_mysolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(CLEmpiricFlux_template
	PUBLIC Eigen3::Eigen
)

target_link_libraries(CLEmpiricFlux_mastersolution.solution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(CLEmpiricFlux_mysolution.solution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(CLEmpiricFlux_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	CLEmpiricFlux_mastersolution.solution
)

target_link_libraries(CLEmpiricFlux_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	CLEmpiricFlux_mysolution.solution
)

gtest_discover_tests(CLEmpiricFlux_test_mysolution)
gtest_discover_tests(CLEmpiricFlux_test_mastersolution)
