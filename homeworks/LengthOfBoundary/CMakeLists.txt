include(GoogleTest)

set(mastersolution_sources
mastersolution/boundarylength.h 
mastersolution/boundarylength.cc 
mastersolution/boundarylength_main.cc)

set(mastersolution_sources
mastersolution/boundarylength.h 
mastersolution/boundarylength.cc 
mastersolution/boundarylength_main.cc)

set(test_master_sources
test/master_tests.cc)

add_library(LengthOfBoundary_master.solution ${mastersolution_sources})
add_executable(LengthOfBoundary_test_mastersolution ${test_master_sources})

add_executable(LengthOfBoundary_mastersolution ${mastersolution_sources})

target_link_libraries(LengthOfBoundary_mastersolution 
    Eigen3::Eigen
    Boost::filesystem
    Boost::system
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
    Boost::filesystem
    Boost::system
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

