include(GoogleTest)

set(template_sources
  templates/1dwaveabsorbingbc_main.cc
  templates/1dwaveabsorbingbc.h
  templates/1dwaveabsorbingbc.cc
)
set(mastersolution_sources
  mastersolution/1dwaveabsorbingbc_main.cc
  mastersolution/1dwaveabsorbingbc.h
  mastersolution/1dwaveabsorbingbc.cc
)
configure_file(mastersolution/viswave.py viswave_mastersolution.py COPYONLY)
configure_file(mastersolution/visenergies.py visenergies_mastersolution.py COPYONLY)

set(mysolution_sources
  mysolution/1dwaveabsorbingbc_main.cc
  mysolution/1dwaveabsorbingbc.h
  mysolution/1dwaveabsorbingbc.cc
)

set(mastersolution_test_sources
  test/1dwaveabsorbingbc_test_mastersolution.cc
)

set(mysolution_test_sources
  test/1dwaveabsorbingbc_test_mysolution.cc
)

# Libraries
add_library(1DWaveAbsorbingBC_mastersolution.solution ${mastersolution_sources})

add_library(1DWaveAbsorbingBC_mysolution.solution ${mysolution_sources})

# Executables
add_executable(1DWaveAbsorbingBC_template ${template_sources})

add_executable(1DWaveAbsorbingBC_mysolution ${mysolution_sources})

add_executable(1DWaveAbsorbingBC_mastersolution ${mastersolution_sources})

add_executable(1DWaveAbsorbingBC_test_mastersolution ${mastersolution_test_sources})

add_executable(1DWaveAbsorbingBC_test_mysolution ${mysolution_test_sources})


target_link_libraries(1DWaveAbsorbingBC_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(1DWaveAbsorbingBC_mysolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(1DWaveAbsorbingBC_template
	PUBLIC Eigen3::Eigen
)

target_link_libraries(1DWaveAbsorbingBC_mastersolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(1DWaveAbsorbingBC_mysolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(1DWaveAbsorbingBC_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	1DWaveAbsorbingBC_mastersolution.solution
)

target_link_libraries(1DWaveAbsorbingBC_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	1DWaveAbsorbingBC_mysolution.solution
)

gtest_discover_tests(1DWaveAbsorbingBC_test_mysolution)
gtest_discover_tests(1DWaveAbsorbingBC_test_mastersolution)
