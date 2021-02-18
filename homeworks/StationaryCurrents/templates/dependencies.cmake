# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/stationarycurrents_main.cc
  ${DIR}/stationarycurrents.h
  ${DIR}/stationarycurrents.cc
  ${DIR}/stationarycurrents_supplement.h
  ${DIR}/stationarycurrents_supplement.cc
)

set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.utils
  LF::lf.uscalfe
)
