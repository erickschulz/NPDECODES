include(GoogleTest)

set(mastersolution_sources
  mastersolution/trace_error_estimate_main.cc
  mastersolution/tee_lapl_robin_assembly.h
  mastersolution/tee_lapl_robin_assembly.cc
)

set(mysolution_sources
  mysolution/trace_error_estimate_main.cc
  mysolution/tee_lapl_robin_assembly.h
  mysolution/tee_lapl_robin_assembly.cc
)

set(template_sources
  templates/trace_error_estimate_main.cc
  templates/tee_lapl_robin_assembly.h
  templates/tee_lapl_robin_assembly.cc
)

set(mastersolution_test_sources
  test/tee_lapl_robin_assembly_test_mastersolution.cc
)


set(mysolution_test_sources
  test/tee_lapl_robin_assembly_test_mysolution.cc
)

# Libraries
add_library(ErrorEstimatesForTraces_mastersolution.solution ${mastersolution_sources})
add_library(ErrorEstimatesForTraces_mysolution.solution ${mysolution_sources})

# Executables
add_executable(ErrorEstimatesForTraces_mastersolution ${mastersolution_sources})
add_executable(ErrorEstimatesForTraces_mysolution ${mysolution_sources})
add_executable(ErrorEstimatesForTraces_template ${template_sources})
add_executable(ErrorEstimatesForTraces_test_mastersolution ${mastersolution_test_sources})
add_executable(ErrorEstimatesForTraces_test_mysolution ${mysolution_test_sources})


target_link_libraries(ErrorEstimatesForTraces_mastersolution
	PUBLIC Eigen3::Eigen
	Boost::system
        Boost::filesystem
        Boost::program_options
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
)

target_link_libraries(ErrorEstimatesForTraces_mysolution
	PUBLIC Eigen3::Eigen
	Boost::system
        Boost::filesystem
        Boost::program_options
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
)

target_link_libraries(ErrorEstimatesForTraces_template
	PUBLIC Eigen3::Eigen
	Boost::system
        Boost::filesystem
        Boost::program_options
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
)
target_link_libraries(ErrorEstimatesForTraces_mastersolution.solution
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
)

target_link_libraries(ErrorEstimatesForTraces_mysolution.solution
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
)

target_link_libraries(ErrorEstimatesForTraces_test_mastersolution
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
        ErrorEstimatesForTraces_mastersolution.solution
)

target_link_libraries(ErrorEstimatesForTraces_test_mysolution
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
        ErrorEstimatesForTraces_mysolution.solution
)


gtest_discover_tests(ErrorEstimatesForTraces_test_mysolution)
gtest_discover_tests(ErrorEstimatesForTraces_test_mastersolution)

