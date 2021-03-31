# Dependencies of mastersolution:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/odesolve_main.cc
  ${DIR}/odesolve.h
  ${DIR}/odesolve.cc
  ${DIR}/polyfit.h
  ${DIR}/polyfit.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
