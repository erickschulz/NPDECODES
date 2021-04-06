# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/exponentialintegrator_main.cc
  ${DIR}/exponentialintegrator.h
  ${DIR}/exponentialintegrator.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
