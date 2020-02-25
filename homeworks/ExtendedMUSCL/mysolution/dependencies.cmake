# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/extendedmuscl.h
  ${DIR}/extendedmuscl.cc
  ${DIR}/slopelimfluxdiff.h
  ${DIR}/sspdriver_main.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
