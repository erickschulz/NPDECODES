# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/systemode_main.cc
  ${DIR}/systemode.h
  ${DIR}/polyfit.h
)

set(LIBRARIES
  Eigen3::Eigen
)
