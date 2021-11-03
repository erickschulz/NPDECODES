# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/fluxlimitedfv_main.cc
  ${DIR}/fluxlimitedfv.h
  ${DIR}/fluxlimitedfv.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
