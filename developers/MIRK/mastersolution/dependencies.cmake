#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/mirk_main.cc
  ${DIR}/mirk.h
)

set(LIBRARIES
  Eigen3::Eigen
)
