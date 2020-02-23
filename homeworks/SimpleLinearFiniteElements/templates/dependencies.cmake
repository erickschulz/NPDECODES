# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/simplelinearfiniteelements_main.cc
  ${DIR}/simplelinearfiniteelements.h
  ${DIR}/simplelinearfiniteelements.cc
  ${DIR}/tria_mesh_2D.h
  ${DIR}/tria_mesh_2D.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
