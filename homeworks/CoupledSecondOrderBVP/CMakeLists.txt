

set(mastersolution_sources
  mastersolution/coupledsecondorderbvp_main.cc
  mastersolution/coupledsecondorderbvp.h
)

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")


# Executables
add_executable(CoupledSecondOrderBVP_mastersolution ${mastersolution_sources})



target_link_libraries(CoupledSecondOrderBVP_mastersolution
	PUBLIC Eigen3::Eigen
	LF::lf.base
	LF::lf.quad
  LF::lf.uscalfe
  LF::lf.mesh
  LF::lf.mesh.utils
  LF::lf.mesh.test_utils
  LF::lf.mesh.hybrid2d
  LF::lf.refinement
  LF::lf.assemble
  LF::lf.io
)
