# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/TestQuadratureRules.cc
  ${DIR}/test_quad_rules.h
  ${DIR}/test_quad_rules.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.base
  LF::lf.quad
)
