# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/finitevolumesineconslaw_main.cc
  ${DIR}/finitevolumesineconslaw.h
  ${DIR}/finitevolumesineconslaw.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
