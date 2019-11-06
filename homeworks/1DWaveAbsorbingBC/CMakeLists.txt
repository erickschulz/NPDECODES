include(GoogleTest)

set(mastersolution_sources
  mastersolution/1dwaveabsorbingbc_main.cc
  mastersolution/1dwaveabsorbingbc.h
  mastersolution/1dwaveabsorbingbc.cc
)
configure_file(mastersolution/viswave.py viswave_mastersolution.py COPYONLY)
configure_file(mastersolution/visenergies.py visenergies_mastersolution.py COPYONLY)

set(mastersolution_test_sources
  test/1dwaveabsorbingbc_test_mastersolution.cc
)

# Libraries
add_library(1DWaveAbsorbingBC_mastersolution.solution ${mastersolution_sources})

# Executables
add_executable(1DWaveAbsorbingBC_mastersolution ${mastersolution_sources})
add_executable(1DWaveAbsorbingBC_test_mastersolution ${mastersolution_test_sources})


# Linking
target_link_libraries(1DWaveAbsorbingBC_mastersolution
	PUBLIC Eigen3::Eigen
)

target_link_libraries(1DWaveAbsorbingBC_mastersolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(1DWaveAbsorbingBC_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::gtest_main
	1DWaveAbsorbingBC_mastersolution.solution
)

gtest_discover_tests(1DWaveAbsorbingBC_test_mastersolution)
