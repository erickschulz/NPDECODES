# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/initcondlv_main.cc
  ${DIR}/initcondlv.h
  ${DIR}/initcondlv.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
