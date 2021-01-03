set(SOURCES
  ${DIR}/expfittedupwind_main.cc
  ${DIR}/expfittedupwind.h
  ${DIR}/expfittedupwind.cc
)

set(LIBRARIES
  Eigen3::Eigen
	LF::lf.assemble
  LF::lf.base
	LF::lf.geometry
	LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.uscalfe
  LF::lf.refinement
)
