#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/ordnotall_main.cc
  ${DIR}/ordnotall.h
  ${DIR}/ordnotall.cc
  ${DIR}/rkintegrator.h
)

set(LIBRARIES
  Eigen3::Eigen
)
