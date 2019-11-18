#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/linearfe1d_main.cc
  ${DIR}/linearfe1d.h
)

set(LIBRARIES
  Eigen3::Eigen
)
