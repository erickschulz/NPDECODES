include(GoogleTest)

set(template_sources
  templates/engquistoshernumericalflux_main.cc
  templates/engquistoshernumericalflux.h
  templates/engquistoshernumericalflux.cc
)
set(mastersolution_sources
  mastersolution/engquistoshernumericalflux_main.cc
  mastersolution/engquistoshernumericalflux.h
  mastersolution/engquistoshernumericalflux.cc
)
configure_file(mastersolution/plot.py plot_mastersolution.py COPYONLY)

set(mysolution_sources
  mysolution/engquistoshernumericalflux_main.cc
  mysolution/engquistoshernumericalflux.h
  mysolution/engquistoshernumericalflux.cc
)

set(mastersolution_test_sources
  test/engquistoshernumericalflux_test_mastersolution.cc
)

set(mysolution_test_sources
  test/engquistoshernumericalflux_test_mysolution.cc
)

# Libraries
add_library(EngquistOsherNumericalFlux_mastersolution.solution ${mastersolution_sources})

add_library(EngquistOsherNumericalFlux_mysolution.solution ${mysolution_sources})

# Executables
add_executable(EngquistOsherNumericalFlux_template ${template_sources})

add_executable(EngquistOsherNumericalFlux_mysolution ${mysolution_sources})

add_executable(EngquistOsherNumericalFlux_mastersolution ${mastersolution_sources})

add_executable(EngquistOsherNumericalFlux_test_mastersolution ${mastersolution_test_sources})

add_executable(EngquistOsherNumericalFlux_test_mysolution ${mysolution_test_sources})


target_link_libraries(EngquistOsherNumericalFlux_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(EngquistOsherNumericalFlux_mysolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(EngquistOsherNumericalFlux_template
	PUBLIC Eigen3::Eigen
)

target_link_libraries(EngquistOsherNumericalFlux_mastersolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(EngquistOsherNumericalFlux_mysolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(EngquistOsherNumericalFlux_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	EngquistOsherNumericalFlux_mastersolution.solution
)

target_link_libraries(EngquistOsherNumericalFlux_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	EngquistOsherNumericalFlux_mysolution.solution
)

gtest_discover_tests(EngquistOsherNumericalFlux_test_mysolution)
gtest_discover_tests(EngquistOsherNumericalFlux_test_mastersolution)
