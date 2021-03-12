#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/initcondlv_main.cc
  ${DIR}/initcondlv.h
  ${DIR}/ode45.h
)

set(LIBRARIES
  Eigen3::Eigen
)
