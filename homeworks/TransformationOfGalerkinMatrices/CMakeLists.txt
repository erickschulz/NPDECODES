
include(GoogleTest)

set(template_sources
  templates/trans_gal_mat.cc
  templates/trans_gal_mat.h
  templates/TransformationOfGalerkinMatrices.cc
)

set(mastersolution_sources
  mastersolution/trans_gal_mat.cc
  mastersolution/trans_gal_mat.h
  mastersolution/TransformationOfGalerkinMatrices.cc
)

set(mysolution_sources
  mysolution/trans_gal_mat.cc
  mysolution/trans_gal_mat.h
  mysolution/TransformationOfGalerkinMatrices.cc
)

set(mastersolution_test_sources
  test/trans_gal_mat_test_mastersolution.cc
)

set(mysolution_test_sources
  test/trans_gal_mat_test_mysolution.cc
)

# Libraries
add_library(TransformationOfGalerkinMatrices_mastersolution.solution ${mastersolution_sources})

add_library(TransformationOfGalerkinMatrices_mysolution.solution ${mysolution_sources})

# Executables
add_executable(TransformationOfGalerkinMatrices_template ${template_sources})

add_executable(TransformationOfGalerkinMatrices_mysolution ${mysolution_sources})

add_executable(TransformationOfGalerkinMatrices_mastersolution ${mastersolution_sources})

add_executable(TransformationOfGalerkinMatrices_test_mastersolution ${mastersolution_test_sources})

add_executable(TransformationOfGalerkinMatrices_test_mysolution ${mysolution_test_sources})


target_link_libraries(TransformationOfGalerkinMatrices_mastersolution
	PUBLIC Eigen3::Eigen
	Boost::system
)

target_link_libraries(TransformationOfGalerkinMatrices_mysolution
	PUBLIC Eigen3::Eigen
	Boost::system
)

target_link_libraries(TransformationOfGalerkinMatrices_template
	PUBLIC Eigen3::Eigen
	Boost::system
)

target_link_libraries(TransformationOfGalerkinMatrices_mastersolution.solution
  PUBLIC Eigen3::Eigen 
  Boost::system
)

target_link_libraries(TransformationOfGalerkinMatrices_mysolution.solution
  PUBLIC Eigen3::Eigen 
  Boost::system
)

target_link_libraries(TransformationOfGalerkinMatrices_test_mastersolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	TransformationOfGalerkinMatrices_mastersolution.solution
	Boost::system
)

target_link_libraries(TransformationOfGalerkinMatrices_test_mysolution 
	PUBLIC Eigen3::Eigen	
	GTest::main
	TransformationOfGalerkinMatrices_mysolution.solution 
	Boost::system

)

gtest_discover_tests(TransformationOfGalerkinMatrices_test_mysolution)
gtest_discover_tests(TransformationOfGalerkinMatrices_test_mastersolution)
