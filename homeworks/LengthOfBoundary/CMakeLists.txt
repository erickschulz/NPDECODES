include(GoogleTest)

set(mastersolution_sources
mastersolution/boundarylength.h 
mastersolution/boundarylength.cc 
mastersolution/boundarylength_main.cc)

set(mysolution_sources
mysolution/boundarylength.h 
mysolution/boundarylength.cc 
mysolution/boundarylength_main.cc)

set(template_sources
templates/boundarylength.h 
templates/boundarylength.cc 
templates/boundarylength_main.cc)

set(test_my_sources
test/tests.cc)
set(test_master_sources
test/master_tests.cc)

add_library(LengthOfBoundary.solution ${mysolution_sources})
add_library(LengthOfBoundary_master.solution ${mastersolution_sources})
add_executable(LengthOfBoundary_test_mysolution ${test_my_sources})
add_executable(LengthOfBoundary_test_mastersolution ${test_master_sources})

add_executable(LengthOfBoundary_mastersolution ${mastersolution_sources})
add_executable(LengthOfBoundary_mysolution ${mysolution_sources})
add_executable(LengthOfBoundary_template ${template_sources})

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

target_link_libraries(LengthOfBoundary_mysolution 
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

target_link_libraries(LengthOfBoundary_template 
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

    target_link_libraries(LengthOfBoundary.solution 
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
    
target_link_libraries(LengthOfBoundary_test_mysolution 
    GTest::main
    LengthOfBoundary.solution
    )
target_link_libraries(LengthOfBoundary_test_mastersolution 
    GTest::main
    LengthOfBoundary_master.solution
    )


add_custom_target(LengthOfBoundary_maketemplate
    COMMAND ${CMAKE_SOURCE_DIR}/remove-lines.sh 
)

gtest_discover_tests(LengthOfBoundary_test_mysolution)
gtest_discover_tests(LengthOfBoundary_test_mastersolution)

