#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/contourplot_main.cc
  ${DIR}/contourplot.h
  ${DIR}/contourplot.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
