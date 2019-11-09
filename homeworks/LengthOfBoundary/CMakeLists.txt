include(GoogleTest)

set(mastersolution_sources
    mastersolution/boundarylength.h 
    mastersolution/boundarylength.cc 
    mastersolution/boundarylength_main.cc
)

set(test_master_sources test/master_tests.cc)

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

add_library(LengthOfBoundary_master.solution ${mastersolution_sources})

add_executable(LengthOfBoundary_mastersolution ${mastersolution_sources})
add_executable(LengthOfBoundary_test_mastersolution ${test_master_sources})

target_link_libraries(LengthOfBoundary_mastersolution 
    Eigen3::Eigen
    LF::lf.base
    LF::lf.mesh
    LF::lf.mesh.utils
    LF::lf.mesh.test_utils
    LF::lf.mesh.hybrid2d
    LF::lf.refinement
    LF::lf.assemble
    LF::lf.io
    LF::lf.uscalfe
)

target_link_libraries(LengthOfBoundary_master.solution 
    Eigen3::Eigen
    LF::lf.base
    LF::lf.mesh
    LF::lf.mesh.utils
    LF::lf.mesh.test_utils
    LF::lf.mesh.hybrid2d
    LF::lf.refinement
    LF::lf.assemble
    LF::lf.io
    LF::lf.uscalfe
)

target_link_libraries(LengthOfBoundary_test_mastersolution 
    GTest::gtest_main
    LengthOfBoundary_master.solution
)

gtest_discover_tests(LengthOfBoundary_test_mastersolution)

