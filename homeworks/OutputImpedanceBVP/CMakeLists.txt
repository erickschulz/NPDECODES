include(GoogleTest)

set(mastersolution_sources
  mastersolution/OutputImpedanceBVP_main.cc
  mastersolution/OutputImpedanceBVP.cc
  mastersolution/OutputImpedanceBVP.h
  mastersolution/evalclass.h
  mastersolution/evalclass.cc
)

set(mysolution_sources
  mysolution/OutputImpedanceBVP_main.cc
  mysolution/OutputImpedanceBVP.cc
  mysolution/OutputImpedanceBVP.h
  mysolution/evalclass.h
  mysolution/evalclass.cc
)

set(template_sources
  templates/OutputImpedanceBVP_main.cc
  templates/OutputImpedanceBVP.cc
  templates/OutputImpedanceBVP.h
  templates/evalclass.h
  templates/evalclass.cc
)

set(mastersolution_test_sources
  test/OutputImpedanceBVP_test_mastersolution.cc
)

set(mysolution_test_sources
  test/OutputImpedanceBVP_test_mysolution.cc
)

# Libraries
add_library(OutputImpedanceBVP_mastersolution.solution ${mastersolution_sources})
add_library(OutputImpedanceBVP_mysolution.solution ${mysolution_sources})

# Executables
add_executable(OutputImpedanceBVP_mastersolution ${mastersolution_sources})
add_executable(OutputImpedanceBVP_mysolution ${mysolution_sources})
add_executable(OutputImpedanceBVP_template ${template_sources})
add_executable(OutputImpedanceBVP_test_mastersolution ${mastersolution_test_sources})
add_executable(OutputImpedanceBVP_test_mysolution ${mysolution_test_sources})


target_link_libraries(OutputImpedanceBVP_mastersolution
	PUBLIC Eigen3::Eigen
	Boost::system
	Boost::filesystem
	LF::lf.base
	LF::lf.uscalfe
	LF::lf.mesh
	LF::lf.mesh.utils
	LF::lf.mesh.hybrid2d
	LF::lf.io
)

target_link_libraries(OutputImpedanceBVP_mysolution
	PUBLIC Eigen3::Eigen
	Boost::system
	Boost::filesystem
	LF::lf.base
	LF::lf.uscalfe
	LF::lf.mesh
	LF::lf.mesh.utils
	LF::lf.mesh.hybrid2d
	LF::lf.io
)

target_link_libraries(OutputImpedanceBVP_template
	PUBLIC Eigen3::Eigen
	Boost::system
	Boost::filesystem
	LF::lf.base
	LF::lf.uscalfe
	LF::lf.mesh
	LF::lf.mesh.utils
	LF::lf.mesh.hybrid2d
	LF::lf.io
)
target_link_libraries(OutputImpedanceBVP_mastersolution.solution
	PUBLIC Eigen3::Eigen
	Boost::system
	Boost::filesystem
	LF::lf.base
	LF::lf.uscalfe
	LF::lf.mesh
	LF::lf.mesh.utils
	LF::lf.mesh.hybrid2d
	LF::lf.io
)

target_link_libraries(OutputImpedanceBVP_mysolution.solution
	PUBLIC Eigen3::Eigen
	Boost::system
	Boost::filesystem
	LF::lf.base
	LF::lf.uscalfe
	LF::lf.mesh
	LF::lf.mesh.utils
	LF::lf.mesh.hybrid2d
	LF::lf.io
)

target_link_libraries(OutputImpedanceBVP_test_mastersolution
	PUBLIC Eigen3::Eigen
	Boost::system
        Boost::filesystem
        Boost::program_options
        GTest::main
	LF::lf.base
	LF::lf.quad
        LF::lf.uscalfe
        LF::lf.mesh
        LF::lf.mesh.utils
        LF::lf.mesh.test_utils
        LF::lf.mesh.hybrid2d
        LF::lf.refinement
        LF::lf.assemble
        LF::lf.io
        OutputImpedanceBVP_mastersolution.solution
)

target_link_libraries(OutputImpedanceBVP_test_mysolution
	PUBLIC Eigen3::Eigen
	Boost::system
        Boost::filesystem
        Boost::program_options
        GTest::main
	LF::lf.base
	LF::lf.quad
        LF::lf.uscalfe
        LF::lf.mesh
        LF::lf.mesh.utils
        LF::lf.mesh.test_utils
        LF::lf.mesh.hybrid2d
        LF::lf.refinement
        LF::lf.assemble
        LF::lf.io
        OutputImpedanceBVP_mysolution.solution
)


gtest_discover_tests(OutputImpedanceBVP_test_mysolution)
gtest_discover_tests(OutputImpedanceBVP_test_mastersolution)

