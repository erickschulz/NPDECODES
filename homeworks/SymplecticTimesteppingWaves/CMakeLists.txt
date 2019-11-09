include(GoogleTest)

set(mastersolution_sources
  mastersolution/symplectictimesteppingwaves_main.cc
  mastersolution/symplectictimesteppingwaves_ode.h
  mastersolution/symplectictimesteppingwaves_ode.cc
  mastersolution/symplectictimesteppingwaves.h
  mastersolution/symplectictimesteppingwaves.cc
  mastersolution/symplectictimesteppingwaves_assemble.h
  )
set(templates_sources
  templates/symplectictimesteppingwaves_main.cc
  templates/symplectictimesteppingwaves_ode.h
  templates/symplectictimesteppingwaves_ode.cc
  templates/symplectictimesteppingwaves.h
  templates/symplectictimesteppingwaves.cc
  templates/symplectictimesteppingwaves_assemble.h
  )
set(mysolution_sources
  mysolution/symplectictimesteppingwaves_main.cc
  mysolution/symplectictimesteppingwaves_ode.h
  mysolution/symplectictimesteppingwaves_ode.cc
  mysolution/symplectictimesteppingwaves.h
  mysolution/symplectictimesteppingwaves.cc
  mysolution/symplectictimesteppingwaves_assemble.h
  )

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")


# Executables
add_executable(symplectictimesteppingwaves_mastersolution ${mastersolution_sources})
add_executable(symplectictimesteppingwaves_templates ${templates_sources})
add_executable(symplectictimesteppingwaves_mysolution ${mysolution_sources})


target_link_libraries(symplectictimesteppingwaves_mastersolution
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

target_link_libraries(symplectictimesteppingwaves_templates
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
target_link_libraries(symplectictimesteppingwaves_mysolution
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
