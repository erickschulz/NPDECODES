#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/exponentialintegrator_main.cc
  ${DIR}/exponentialintegrator.h
  ${DIR}/exponentialintegrator.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
