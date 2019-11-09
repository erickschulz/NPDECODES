include(GoogleTest)

set(mastersolution_sources
  mastersolution/BoundaryWave_main.cc
  mastersolution/BoundaryWave.cc
  mastersolution/BoundaryWave.h
)

 set(mysolution_sources
   mysolution/BoundaryWave_main.cc
   mysolution/BoundaryWave.cc
   mysolution/BoundaryWave.h
 )

 set(template_sources
   templates/BoundaryWave_main.cc
   templates/BoundaryWave.cc
   templates/BoundaryWave.h
 )

set(mastersolution_test_sources
  test/BoundaryWave_test_mastersolution.cc
)

set(mastersolution_test_sources
  test/BoundaryWave_test_mastersolution.cc
)

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Libraries
add_library(BoundaryWave_mastersolution.solution ${mastersolution_sources})
add_library(BoundaryWave_mysolution.solution ${mysolution_sources})

# Executables
add_executable(BoundaryWave_mastersolution ${mastersolution_sources})
 add_executable(BoundaryWave_mysolution ${mysolution_sources})
 add_executable(BoundaryWave_template ${template_sources})
#add_executable(BoundaryWave_test_mastersolution ${mastersolution_test_sources})
#add_executable(BoundaryWave_test_mysolution ${mysolution_test_sources})


target_link_libraries(BoundaryWave_mastersolution
	PUBLIC Eigen3::Eigen
	LF::lf.base
	LF::lf.uscalfe
	LF::lf.mesh
	LF::lf.mesh.utils
	LF::lf.mesh.hybrid2d
	LF::lf.io
)

target_link_libraries(BoundaryWave_mysolution
 	PUBLIC Eigen3::Eigen
 	LF::lf.base
 	LF::lf.uscalfe
 	LF::lf.mesh
 	LF::lf.mesh.utils
 	LF::lf.mesh.hybrid2d
 	LF::lf.io
 )

target_link_libraries(BoundaryWave_template
 	PUBLIC Eigen3::Eigen
 	LF::lf.base
 	LF::lf.uscalfe
 	LF::lf.mesh
 	LF::lf.mesh.utils
 	LF::lf.mesh.hybrid2d
 	LF::lf.io
)

target_link_libraries(BoundaryWave_mastersolution.solution
	PUBLIC Eigen3::Eigen
	LF::lf.base
	LF::lf.uscalfe
	LF::lf.mesh
	LF::lf.mesh.utils
	LF::lf.mesh.hybrid2d
	LF::lf.io
)

 target_link_libraries(BoundaryWave_mysolution.solution
 	PUBLIC Eigen3::Eigen
 	LF::lf.base
 	LF::lf.uscalfe
 	LF::lf.mesh
 	LF::lf.mesh.utils
 	LF::lf.mesh.hybrid2d
 	LF::lf.io
)

#target_link_libraries(BoundaryWave_test_mastersolution
#	PUBLIC Eigen3::Eigen
#        GTest::main
#	LF::lf.base
#	LF::lf.quad
#        LF::lf.uscalfe
#        LF::lf.mesh
#        LF::lf.mesh.utils
#        LF::lf.mesh.test_utils
#        LF::lf.mesh.hybrid2d
#        LF::lf.refinement
#        LF::lf.assemble
#        LF::lf.io
#        BoundaryWave_mastersolution.solution
#)

#target_link_libraries(BoundaryWave_test_mysolution
#	PUBLIC Eigen3::Eigen
#        GTest::main
#	LF::lf.base
#	LF::lf.quad
#        LF::lf.uscalfe
#        LF::lf.mesh
#        LF::lf.mesh.utils
#        LF::lf.mesh.test_utils
#        LF::lf.mesh.hybrid2d
#        LF::lf.refinement
#        LF::lf.assemble
#        LF::lf.io
#        BoundaryWave_mysolution.solution
#)

#gtest_discover_tests(BoundaryWave_test_mysolution)
#gtest_discover_tests(BoundaryWave_test_mastersolution)

