# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/engquistoshernumericalflux_main.cc
  ${DIR}/engquistoshernumericalflux.h
  ${DIR}/engquistoshernumericalflux.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
