# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/parametricelementmatrices_main.cc
  ${DIR}/anisotropicdiffusionelementmatrixprovider.h
  ${DIR}/anisotropicdiffusionelementmatrixprovider.cc
  ${DIR}/fesourceelemvecprovider.h
  ${DIR}/fesourceelemvecprovider.cc
  ${DIR}/impedanceboundaryedgematrixprovider.h
  ${DIR}/impedanceboundaryedgematrixprovider.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.geometry
  LF::lf.io
  LF::lf.mesh
  LF::lf.mesh.test_utils
  LF::lf.quad
  LF::lf.uscalfe
)
