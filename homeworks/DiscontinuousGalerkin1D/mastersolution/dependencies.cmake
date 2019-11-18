# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/discontinuousgalerkin1d_main.cc
  ${DIR}/discontinuousgalerkin1d.h
  ${DIR}/discontinuousgalerkin1d.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
