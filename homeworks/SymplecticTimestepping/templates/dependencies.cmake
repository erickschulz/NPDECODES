# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/symplectictimestepping_main.cc
  ${DIR}/symplectictimestepping.h
  ${DIR}/symplectictimestepping.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
