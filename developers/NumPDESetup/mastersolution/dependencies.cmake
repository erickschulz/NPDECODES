#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/numpdesetup_main.cc
  ${DIR}/numpdesetup.h
  ${DIR}/numpdesetup.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
