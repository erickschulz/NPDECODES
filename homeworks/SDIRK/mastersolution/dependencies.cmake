# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/sdirk_main.cc
  ${DIR}/sdirk.h
  ${DIR}/sdirk.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
