#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/odesolve_main.cc
  ${DIR}/odesolve.h
)

set(LIBRARIES
  Eigen3::Eigen
)
