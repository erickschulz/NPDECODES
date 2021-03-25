# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/contourplot_main.cc
  ${DIR}/contourplot.h
  ${DIR}/contourplot.cc
  ${DIR}/ode45.h
)

set(LIBRARIES
  Eigen3::Eigen
)
