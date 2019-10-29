include(GoogleTest)

set(libraries
  Eigen3::Eigen
  Boost::filesystem 
  Boost::system 
  Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.geometry
  LF::lf.io
  LF::lf.quad
  LF::lf.uscalfe
)

set(template_sources
  template/ElementMatrixComputation_main.cc
  template/MyLinearFEElementMatrix.h
  template/MyLinearFEElementMatrix.cc
  template/MyLinearLoadVector.h
  template/MyLinearLoadVector.cc
  template/solve.h
  meshes/mesh.h
)

set(mysolution_sources
  mysolution/MyLinearFEElementMatrix.h
  mysolution/MyLinearFEElementMatrix.cc
  mysolution/MyLinearLoadVector.h
  mysolution/MyLinearLoadVector.cc
  mysolution/solve.h
  meshes/mesh.h
)

set(mastersolution_sources
  mastersolution/MyLinearFEElementMatrix.h
  mastersolution/MyLinearFEElementMatrix.cc
  mastersolution/MyLinearLoadVector.h
  mastersolution/MyLinearLoadVector.cc
  mastersolution/solve.h
  meshes/mesh.h
)

set(test_sources
  test/ElementMatrixComputation_tests.cc
)

add_library(ElementMatrixComputation.mysolution ${mysolution_sources})
target_include_directories(ElementMatrixComputation.mysolution PUBLIC ./mysolution)
target_link_libraries(ElementMatrixComputation.mysolution PUBLIC ${libraries})

add_library(ElementMatrixComputation.mastersolution ${mastersolution_sources})
# specify the include directory of the master solution library
# The PUBLIC specifier is used so that dependent executables/libraries also 
# use this include directory.
target_include_directories(ElementMatrixComputation.mastersolution PUBLIC ./mastersolution)

# Link the mastersolution library with all LehrFEM++/boost libraries$
# The PUBLIC specifier is used so that all dependent executables/libraries that
# use this library also link with the `${libraries}`
target_link_libraries(ElementMatrixComputation.mastersolution PUBLIC ${libraries})

add_executable(ElementMatrixComputation_template ${template_sources})
add_executable(ElementMatrixComputation_mysolution ${mysolution_sources} mysolution/ElementMatrixComputation_main.cc)
add_executable(ElementMatrixComputation_mastersolution ${mastersolution_sources} mastersolution/ElementMatrixComputation_main.cc)
add_executable(ElementMatrixComputation_test_mysolution ${test_sources})
add_executable(ElementMatrixComputation_test_mastersolution ${test_sources})

target_link_libraries(ElementMatrixComputation_template PUBLIC ${libraries})
target_link_libraries(ElementMatrixComputation_mastersolution ElementMatrixComputation.mastersolution)
target_link_libraries(ElementMatrixComputation_mysolution ElementMatrixComputation.mysolution)

target_link_libraries(ElementMatrixComputation_test_mysolution GTest::main ElementMatrixComputation.mysolution)
target_link_libraries(ElementMatrixComputation_test_mastersolution GTest::main ElementMatrixComputation.mastersolution)

gtest_discover_tests(ElementMatrixComputation_test_mysolution)
gtest_discover_tests(ElementMatrixComputation_test_mastersolution)

add_custom_target(ElementMatrixComputation_maketemplate
    COMMAND ${CMAKE_SOURCE_DIR}/strip-comments.sh
)
