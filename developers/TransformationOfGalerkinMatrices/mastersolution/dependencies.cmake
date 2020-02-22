#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/transformationofgalerkinmatrices_main.cc
  ${DIR}/transformationofgalerkinmatrices.h
  ${DIR}/transformationofgalerkinmatrices.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
