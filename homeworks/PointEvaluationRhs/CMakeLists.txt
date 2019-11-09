include(GoogleTest)


set(mastersolution_sources
mastersolution/norms.cc
mastersolution/pointEvaluation.cc
mastersolution/PointEvaluationRhs_main.cc
)

add_executable(PointEvaluationRhs_mastersolution ${mastersolution_sources})

target_link_libraries(PointEvaluationRhs_mastersolution
  PUBLIC Eigen3::Eigen
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.assemble
  LF::lf.refinement
  LF::lf.io
  LF::lf.uscalfe
  )
