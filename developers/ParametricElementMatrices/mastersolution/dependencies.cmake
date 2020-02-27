#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/parametricelementmatrices_main.cc
  ${DIR}/ansiotropicdiffusionelementmatrixprovider.h
  ${DIR}/ansiotropicdiffusionelementmatrixprovider.cc
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
