#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/fluxlimitedfv_main.cc
  ${DIR}/fluxlimitedfv.h
  ${DIR}/fluxlimitedfv.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
