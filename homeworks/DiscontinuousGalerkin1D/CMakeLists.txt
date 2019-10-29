include(GoogleTest)

set(mastersolution_sources
  mastersolution/discontinuousgalerkin1d_main.cc
  mastersolution/discontinuousgalerkin1d.h
  mastersolution/discontinuousgalerkin1d.cc
)
configure_file(mastersolution/plot.py plot_mastersolution.py COPYONLY)

set(mysolution_sources
  mysolution/discontinuousgalerkin1d_main.cc
  mysolution/discontinuousgalerkin1d.h
  mysolution/discontinuousgalerkin1d.cc
)

set(template_sources
  templates/discontinuousgalerkin1d_main.cc
  templates/discontinuousgalerkin1d.h
  templates/discontinuousgalerkin1d.cc
)

set(mastersolution_test_sources
  test/discontinuousgalerkin1d_test_mastersolution.cc
)

set(mysolution_test_sources
  test/discontinuousgalerkin1d_test_mysolution.cc
)

# Libraries
add_library(DiscontinuousGalerkin1D_mastersolution.solution ${mastersolution_sources})

add_library(DiscontinuousGalerkin1D_mysolution.solution ${mysolution_sources})

# Executables
add_executable(DiscontinuousGalerkin1D_mastersolution ${mastersolution_sources})

add_executable(DiscontinuousGalerkin1D_mysolution ${mysolution_sources})

add_executable(DiscontinuousGalerkin1D_template ${template_sources})

add_executable(DiscontinuousGalerkin1D_test_mastersolution ${mastersolution_test_sources})

add_executable(DiscontinuousGalerkin1D_test_mysolution ${mysolution_test_sources})


target_link_libraries(DiscontinuousGalerkin1D_mastersolution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(DiscontinuousGalerkin1D_mysolution
  PUBLIC Eigen3::Eigen
)

target_link_libraries(DiscontinuousGalerkin1D_template
  PUBLIC Eigen3::Eigen
)

target_link_libraries(DiscontinuousGalerkin1D_mastersolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(DiscontinuousGalerkin1D_mysolution.solution
  PUBLIC Eigen3::Eigen 
)

target_link_libraries(DiscontinuousGalerkin1D_test_mastersolution 
  PUBLIC Eigen3::Eigen	
  GTest::main
  DiscontinuousGalerkin1D_mastersolution.solution
)

target_link_libraries(DiscontinuousGalerkin1D_test_mysolution 
  PUBLIC Eigen3::Eigen	
  GTest::main
  DiscontinuousGalerkin1D_mysolution.solution
)

gtest_discover_tests(DiscontinuousGalerkin1D_test_mysolution)
gtest_discover_tests(DiscontinuousGalerkin1D_test_mastersolution)
