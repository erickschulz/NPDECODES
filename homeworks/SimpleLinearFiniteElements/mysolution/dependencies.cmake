# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/simple_linear_finite_elements_main.cc
  ${DIR}/simple_linear_finite_elements.h
  ${DIR}/simple_linear_finite_elements.cc
  ${DIR}/tria_mesh_2D.cc
  ${DIR}/local_computations.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
