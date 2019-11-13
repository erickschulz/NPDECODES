#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/1dwaveabsorbingbc_main.cc
  ${DIR}/1dwaveabsorbingbc.h
  ${DIR}/1dwaveabsorbingbc.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
