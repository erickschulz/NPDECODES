# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/initcondlv_main.cc
  ${DIR}/initcondlv.h
  ${DIR}/initcondlv.cc
  ${DIR}/ode45.h
)

set(LIBRARIES
  Eigen3::Eigen
)
