# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/nlmatode_main.cc
  ${DIR}/nlmatode.h
  ${DIR}/nlmatode.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
