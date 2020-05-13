# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/extendedmuscl.h
  ${DIR}/extendedmuscl.cc
  ${DIR}/slopelimfluxdiff.h
  ${DIR}/extendedmuscl_main.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
