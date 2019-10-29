
include(GoogleTest)

set(sources
  assembly_tests.cc
  element_tests.cc
  ../SimpleLinearFEM2D.h
  ../local_assembler.h
  ../local_assembler.cc
  ../matrix_assembler.cc
  ../vector_assembler.cc
  ../Mesh.cc
  ../fe_solver.cc
)

add_executable(SimpleLinearFEM2D.test ${sources})
target_link_libraries(SimpleLinearFEM2D.test PUBLIC Eigen3::Eigen
Boost::filesystem Boost::system Boost::program_options 
GTest::main)
gtest_discover_tests(SimpleLinearFEM2D.test)
