#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/taylorode_main.cc
  ${DIR}/taylorode.h
)

set(LIBRARIES
  Eigen3::Eigen
)
