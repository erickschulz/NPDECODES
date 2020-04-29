# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/transformationofgalerkinmatrices_main.cc
  ${DIR}/transformationofgalerkinmatrices.h
  ${DIR}/transformationofgalerkinmatrices.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
