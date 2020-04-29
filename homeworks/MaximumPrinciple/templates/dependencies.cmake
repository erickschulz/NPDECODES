# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/maximumprinciple_main.cc
  ${DIR}/maximumprinciple.h
  ${DIR}/maximumprinciple.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
