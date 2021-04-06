# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/stabrk3_main.cc
  ${DIR}/stabrk3.h
  ${DIR}/stabrk3.cc
  ${DIR}/rkintegrator.h
)

set(LIBRARIES
  Eigen3::Eigen
)
