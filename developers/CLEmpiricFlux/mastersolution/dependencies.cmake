#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/clempiricflux_main.cc
  ${DIR}/clempiricflux.h
  ${DIR}/clempiricflux.cc
  ${DIR}/solvecauchyproblem.h
  ${DIR}/solvecauchyproblem.cc
  ${DIR}/uniformcubicspline.h
  ${DIR}/uniformcubicspline.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
