
include(GoogleTest)

set(template_sources
  template/local_laplace_qfe.h
  template/local_laplace_qfe.cc
  template/main.cc
  template/qfe_interpolator.h
  template/qfe_interpolator.cc
  template/qfe_provider_tester.h
  utils/figure.h
  utils/figure.cc
  utils/MglLabel.h
  utils/MglPlot.h
  utils/MglStyle.h
  utils/MglStyle.h
  utils/FigureConfig.h
)

set(mysolution_sources
  mysolution/local_laplace_qfe.h
  mysolution/local_laplace_qfe.cc
  mysolution/main.cc
  mysolution/qfe_interpolator.h
  mysolution/qfe_interpolator.cc
  mysolution/qfe_provider_tester.h
  utils/figure.h
  utils/figure.cc
  utils/MglLabel.h
  utils/MglPlot.h
  utils/MglStyle.h
  utils/MglStyle.h
  utils/FigureConfig.h
)

set(mastersolution_sources
  mastersolution/local_laplace_qfe.h
  mastersolution/local_laplace_qfe.cc
  mastersolution/main.cc
  mastersolution/qfe_interpolator.h
  mastersolution/qfe_interpolator.cc
  mastersolution/qfe_provider_tester.h
  utils/figure.h
  utils/figure.cc
  utils/MglLabel.h
  utils/MglPlot.h
  utils/MglStyle.h
  utils/MglStyle.h
  utils/FigureConfig.h
)

set(mysolution_test_sources
  test/debugging_fem_test_mysolution.cc
)

set(mastersolution_test_sources
  test/debugging_fem_test_mastersolution.cc
)

# Libraries
add_library(DebuggingFEM_mysolution.solution ${mysolution_sources})

add_library(DebuggingFEM_mastersolution.solution ${mastersolution_sources})

# Executables
add_executable(DebuggingFEM_template ${template_sources})

add_executable(DebuggingFEM_mysolution ${mysolution_sources})

add_executable(DebuggingFEM_test_mysolution ${mysolution_test_sources})

add_executable(DebuggingFEM_mastersolution ${mastersolution_sources})

add_executable(DebuggingFEM_test_mastersolution ${mastersolution_test_sources})

target_link_libraries(DebuggingFEM_template
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  MathGL
)

target_link_libraries(DebuggingFEM_mysolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  MathGL
)

target_link_libraries(DebuggingFEM_mysolution.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  MathGL
)

target_link_libraries(DebuggingFEM_test_mysolution 
  PUBLIC Eigen3::Eigen
  LF::lf.mesh.test_utils
  LF::lf.assemble
  GTest::main
  DebuggingFEM_mysolution.solution
)

target_link_libraries(DebuggingFEM_mastersolution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  MathGL
)

target_link_libraries(DebuggingFEM_mastersolution.solution
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
  LF::lf.uscalfe
  MathGL
)

target_link_libraries(DebuggingFEM_test_mastersolution 
  PUBLIC Eigen3::Eigen
  LF::lf.mesh.test_utils
  LF::lf.assemble
  GTest::main
  DebuggingFEM_mastersolution.solution
)

gtest_discover_tests(DebuggingFEM_test_mysolution)

gtest_discover_tests(DebuggingFEM_test_mastersolution)
