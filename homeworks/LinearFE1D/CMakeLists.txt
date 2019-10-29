include(GoogleTest)

set(template_sources
templates/solve_LinearFE1D.h
templates/LinearFE1D.cc
)

set(mysolution_sources
mysolution/solve_LinearFE1D.h
mysolution/LinearFE1D.cc
)

set(mastersolution_sources
mastersolution/solve_LinearFE1D.h
mastersolution/LinearFE1D.cc
)

set(my_test_sources
test/tests.cc
)
set(master_test_sources
test/master_tests.cc
)


add_executable(LinearFE1D_template ${template_sources})
add_executable(LinearFE1D_mysolution ${mysolution_sources})
add_executable(LinearFE1D_mastersolution ${mastersolution_sources})


add_library(LinearFE1Dmy.solution ${mysolution_sources})
add_library(LinearFE1Dmaster.solution ${mastersolution_sources})
add_executable(LinearFE1D_test_mysolution ${my_test_sources})
add_executable(LinearFE1D_test_mastersolution ${master_test_sources})


target_link_libraries(LinearFE1D_template
Eigen3::Eigen
Boost::filesystem
Boost::system
Boost::program_options
LF::lf.mesh
LF::lf.mesh.utils
LF::lf.mesh.test_utils
LF::lf.mesh.hybrid2d
LF::lf.assemble
LF::lf.io
LF::lf.uscalfe
)

target_link_libraries(LinearFE1D_mysolution
Eigen3::Eigen
Boost::filesystem
Boost::system
Boost::program_options
LF::lf.mesh
LF::lf.mesh.utils
LF::lf.mesh.test_utils
LF::lf.mesh.hybrid2d
LF::lf.assemble
LF::lf.io
LF::lf.uscalfe
)
target_link_libraries(LinearFE1D_mastersolution
Eigen3::Eigen
Boost::filesystem
Boost::system
Boost::program_options
LF::lf.mesh
LF::lf.mesh.utils
LF::lf.mesh.test_utils
LF::lf.mesh.hybrid2d
LF::lf.assemble
LF::lf.io
LF::lf.uscalfe
)

target_link_libraries(LinearFE1Dmy.solution
Eigen3::Eigen
Boost::filesystem
Boost::system
Boost::program_options
LF::lf.mesh
LF::lf.mesh.utils
LF::lf.mesh.test_utils
LF::lf.mesh.hybrid2d
LF::lf.assemble
LF::lf.io
LF::lf.uscalfe
)

target_link_libraries(LinearFE1Dmaster.solution
Eigen3::Eigen
Boost::filesystem
Boost::system
Boost::program_options
LF::lf.mesh
LF::lf.mesh.utils
LF::lf.mesh.test_utils
LF::lf.mesh.hybrid2d
LF::lf.assemble
LF::lf.io
LF::lf.uscalfe
)

target_link_libraries(LinearFE1D_test_mysolution GTest::main LinearFE1Dmy.solution)
target_link_libraries(LinearFE1D_test_mastersolution GTest::main LinearFE1Dmaster.solution)

add_custom_target(LinearFE1D_maketemplate
COMMAND ${CMAKE_SOURCE_DIR}/strip-comments.sh
)

gtest_discover_tests(LinearFE1D_test_mysolution)
gtest_discover_tests(LinearFE1D_test_mastersolution)
