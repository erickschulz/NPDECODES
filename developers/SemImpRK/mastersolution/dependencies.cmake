#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/semimprk_main.cc
  ${DIR}/semimprk.h
  ${DIR}/semimprk.cc
  ${DIR}/polyfit.h
)

set(LIBRARIES
  Eigen3::Eigen
)
