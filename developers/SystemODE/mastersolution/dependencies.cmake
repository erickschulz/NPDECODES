#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/systemode_main.cc
  ${DIR}/systemode.h
  ${DIR}/polyfit.h
  ${DIR}/rkintegrator.h
)

set(LIBRARIES
  Eigen3::Eigen
)
