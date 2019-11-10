include(GoogleTest)

set(mastersolution_sources
  mastersolution/sdirkmethodoflines_main.cc
  mastersolution/sdirkmethodoflines.h
  mastersolution/sdirkmethodoflines.cc
  mastersolution/sdirkmethodoflines_ode.h
  mastersolution/sdirkmethodoflines_ode.cc
)
set(mysolution_sources
  mysolution/sdirkmethodoflines_main.cc
  mysolution/sdirkmethodoflines.h
  mysolution/sdirkmethodoflines.cc
  mysolution/sdirkmethodoflines_ode.h
  mysolution/sdirkmethodoflines_ode.cc
)
set(templates_sources
  templates/sdirkmethodoflines_main.cc
  templates/sdirkmethodoflines.h
  templates/sdirkmethodoflines.cc
  templates/sdirkmethodoflines_ode.h
  templates/sdirkmethodoflines_ode.cc
)

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Executables
add_executable(sdirkmethodoflines_mastersolution ${mastersolution_sources})
add_executable(sdirkmethodoflines_mysolution ${mysolution_sources})
add_executable(sdirkmethodoflines_templates ${templates_sources})

target_link_libraries(sdirkmethodoflines_mastersolution
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
target_link_libraries(sdirkmethodoflines_mysolution
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
target_link_libraries(sdirkmethodoflines_templates
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
