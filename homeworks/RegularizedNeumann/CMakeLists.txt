include(GoogleTest)

set(template_sources
templates/regNeumann.h
templates/RegularizedNeumann_main.cc
)
set(mysolution_sources
mysolution/regNeumann.h
mysolution/RegularizedNeumann_main.cc
)
set(mastersolution_sources
mastersolution/regNeumann.h
mastersolution/RegularizedNeumann_main.cc
)
set(master_test_sources
  tests/tests_master.cc
)
set(my_test_sources
  tests/tests_mysolution.cc
)

# add_executable(RegularizedNeumann_template ${template_sources})
# add_executable(RegularizedNeumann_mysolution ${mysolution_sources})
add_executable(RegularizedNeumann_mastersolution ${mastersolution_sources})

# target_link_libraries(RegularizedNeumann_template
#   PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
#   LF::lf.mesh
#   LF::lf.mesh.utils
#   LF::lf.mesh.test_utils
#   LF::lf.mesh.hybrid2d
#   LF::lf.assemble
#   LF::lf.io
#   LF::lf.uscalfe
# )
# target_link_libraries(RegularizedNeumann_mysolution
#   PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
#   LF::lf.mesh
#   LF::lf.mesh.utils
#   LF::lf.mesh.test_utils
#   LF::lf.mesh.hybrid2d
#   LF::lf.io
#   LF::lf.uscalfe
# )

target_link_libraries(RegularizedNeumann_mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
)

#Unit tests
# add_executable(RegularizedNeumann_test_mysolution ${my_test_sources})
add_executable(RegularizedNeumann_test_mastersolution ${master_test_sources})

# add_library(RegularizedNeumann_my.solution ${mysolution_sources})
add_library(RegularizedNeumann_master.solution ${mastersolution_sources})

# target_link_libraries(RegularizedNeumann_my.solution
#   PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
#   LF::lf.mesh
#   LF::lf.mesh.utils
#   LF::lf.mesh.test_utils
#   LF::lf.mesh.hybrid2d
#   LF::lf.assemble
#   LF::lf.io
#   LF::lf.uscalfe
#   )
target_link_libraries(RegularizedNeumann_master.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  )

# target_link_libraries(RegularizedNeumann_test_mysolution
# GTest::main 
# RegularizedNeumann_my.solution
# )

target_link_libraries(RegularizedNeumann_test_mastersolution
GTest::main 
RegularizedNeumann_master.solution
)

gtest_discover_tests(RegularizedNeumann_test_mastersolution)
# gtest_discover_tests(RegularizedNeumann_test_mysolution)


