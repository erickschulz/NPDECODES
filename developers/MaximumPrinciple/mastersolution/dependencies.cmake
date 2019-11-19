#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/maximum_principle_main.cc
  ${DIR}/maximum_principle.h
  ${DIR}/maximum_principle.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
