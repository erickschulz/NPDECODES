include(GoogleTest)

set(mastersolution_sources
  mastersolution/radauthreetimestepping_main.cc
  mastersolution/radau_three_timestepping_ode.h
  mastersolution/radau_three_timestepping_ode.cc
  mastersolution/radau_three_timestepping.h
  mastersolution/radau_three_timestepping.cc
)



# Executables
add_executable(radauthreetimestepping_mastersolution ${mastersolution_sources})



target_link_libraries(radauthreetimestepping_mastersolution
	PUBLIC Eigen3::Eigen
	Boost::system
        Boost::filesystem
        Boost::program_options
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
