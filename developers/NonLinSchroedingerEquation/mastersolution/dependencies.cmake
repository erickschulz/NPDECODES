#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/nonlinschroedingerequation_main.cc
  ${DIR}/nonlinschroedingerequation.h
  ${DIR}/nonlinschroedingerequation.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
