# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/1dwaveabsorbingbc_main.cc
  ${DIR}/1dwaveabsorbingbc.h
  ${DIR}/1dwaveabsorbingbc.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
