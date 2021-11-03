# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/matode_main.cc
  ${DIR}/matode.h
  ${DIR}/matode.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
