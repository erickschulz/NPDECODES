# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/numpdesetup_main.cc
  ${DIR}/numpdesetup.h
  ${DIR}/numpdesetup.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
