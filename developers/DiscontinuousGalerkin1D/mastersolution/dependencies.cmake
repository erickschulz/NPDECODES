#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/discontinuousgalerkin1d_main.cc
  ${DIR}/discontinuousgalerkin1d.h
  ${DIR}/discontinuousgalerkin1d.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
