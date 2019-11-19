# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/TransformationOfGalerkinMatrices.cc
  ${DIR}/trans_gal_mat.h
  ${DIR}/trans_gal_mat.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
