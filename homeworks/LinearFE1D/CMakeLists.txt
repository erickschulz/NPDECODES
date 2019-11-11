include(GoogleTest)

set(mastersolution_sources
mastersolution/linearfe1d.h
mastersolution/linearfe1d_main.cc
)

set(master_test_sources
test/linearfe1d_test_mastersolution.cc
)


add_executable(LinearFE1D_mastersolution ${mastersolution_sources})


add_library(LinearFE1Dmaster.solution ${mastersolution_sources})
add_executable(LinearFE1D_test_mastersolution ${master_test_sources})


target_link_libraries(LinearFE1D_mastersolution
Eigen3::Eigen
)

target_link_libraries(LinearFE1Dmaster.solution
Eigen3::Eigen
)

target_link_libraries(LinearFE1D_test_mastersolution GTest::gtest_main LinearFE1Dmaster.solution)

add_custom_target(LinearFE1D_maketemplate
COMMAND ${CMAKE_SOURCE_DIR}/strip-comments.sh
)

gtest_discover_tests(LinearFE1D_test_mastersolution)
