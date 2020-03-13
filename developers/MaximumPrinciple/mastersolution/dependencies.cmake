#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/maximumprinciple_main.cc
  ${DIR}/maximumprinciple.h
  ${DIR}/maximumprinciple.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
