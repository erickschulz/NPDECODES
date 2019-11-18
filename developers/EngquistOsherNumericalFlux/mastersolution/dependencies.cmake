#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/engquistoshernumericalflux_main.cc
  ${DIR}/engquistoshernumericalflux.h
  ${DIR}/engquistoshernumericalflux.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
