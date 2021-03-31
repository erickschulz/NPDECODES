# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/crossprod_main.cc
  ${DIR}/crossprod.h
  ${DIR}/crossprod.cc
  ${DIR}/dampnewton.h
  ${DIR}/implicitrkintegrator.h
)

set(LIBRARIES
  Eigen3::Eigen
)
