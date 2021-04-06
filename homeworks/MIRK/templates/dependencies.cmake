# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/mirk_main.cc
  ${DIR}/mirk.h
)

set(LIBRARIES
  Eigen3::Eigen
)
