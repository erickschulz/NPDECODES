# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/taylorode_main.cc
  ${DIR}/taylorode.h
  ${DIR}/taylorode.cc
)


set(LIBRARIES
  Eigen3::Eigen
)
