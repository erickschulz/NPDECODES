#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/burgersequation_main.cc
  ${DIR}/burgersequation.h
  ${DIR}/burgersequation.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
