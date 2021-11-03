#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/rk3prey_main.cc
  ${DIR}/rk3prey.h
)

set(LIBRARIES
  Eigen3::Eigen
)
