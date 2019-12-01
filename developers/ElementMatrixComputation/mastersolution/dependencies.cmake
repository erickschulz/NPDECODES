#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/ElementMatrixComputation_main.cc
  ${DIR}/MyLinearFEElementMatrix.h
  ${DIR}/MyLinearFEElementMatrix.cc
  ${DIR}/MyLinearLoadVector.h
  ${DIR}/MyLinearLoadVector.cc
  ${DIR}/solve.h
  ${DIR}/solve.cc
  meshes/mesh.h
  meshes/mesh.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.uscalfe
)
