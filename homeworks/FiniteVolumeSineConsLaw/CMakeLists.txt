include(GoogleTest)

set(template_sources
  templates/finitevolumesineconslaw_main.cc
  templates/finitevolumesineconslaw.h
  templates/finitevolumesineconslaw.cc
)
set(mastersolution_sources
  mastersolution/finitevolumesineconslaw_main.cc
  mastersolution/finitevolumesineconslaw.h
  mastersolution/finitevolumesineconslaw.cc
)
configure_file(mastersolution/plot.py plot_mastersolution.py COPYONLY)

set(mysolution_sources
  mysolution/finitevolumesineconslaw_main.cc
  mysolution/finitevolumesineconslaw.h
  mysolution/finitevolumesineconslaw.cc
)

set(mastersolution_test_sources
  test/finitevolumesineconslaw_test_mastersolution.cc
)

set(mysolution_test_sources
  test/finitevolumesineconslaw_test_mysolution.cc
)

# Libraries
add_library(FiniteVolumeSineConsLaw_mastersolution.solution ${mastersolution_sources})

add_library(FiniteVolumeSineConsLaw_mysolution.solution ${mysolution_sources})

# Executables
add_executable(FiniteVolumeSineConsLaw_template ${template_sources})

add_executable(FiniteVolumeSineConsLaw_mysolution ${mysolution_sources})

add_executable(FiniteVolumeSineConsLaw_mastersolution ${mastersolution_sources})

add_executable(FiniteVolumeSineConsLaw_test_mastersolution ${mastersolution_test_sources})

add_executable(FiniteVolumeSineConsLaw_test_mysolution ${mysolution_test_sources})


target_link_libraries(FiniteVolumeSineConsLaw_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(FiniteVolumeSineConsLaw_mysolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(FiniteVolumeSineConsLaw_template
	PUBLIC Eigen3::Eigen
)

target_link_libraries(FiniteVolumeSineConsLaw_mastersolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(FiniteVolumeSineConsLaw_mysolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(FiniteVolumeSineConsLaw_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	FiniteVolumeSineConsLaw_mastersolution.solution
)

target_link_libraries(FiniteVolumeSineConsLaw_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	FiniteVolumeSineConsLaw_mysolution.solution
)

gtest_discover_tests(FiniteVolumeSineConsLaw_test_mysolution)
gtest_discover_tests(FiniteVolumeSineConsLaw_test_mastersolution)
