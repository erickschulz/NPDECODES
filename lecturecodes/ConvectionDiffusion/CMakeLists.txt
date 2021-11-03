add_executable(lecturecodes.convectiondiffusion.layer layer_main.cc cd_tools.cc ../../homeworks/UpwindQuadrature/mastersolution/upwindquadrature.cc)
target_compile_definitions(lecturecodes.convectiondiffusion.layer PRIVATE CURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/${DIR}\")
target_compile_definitions(lecturecodes.convectiondiffusion.layer PRIVATE CURRENT_BINARY_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
target_link_libraries(lecturecodes.convectiondiffusion.layer
    PUBLIC Eigen3::Eigen
    LF::lf.assemble
    LF::lf.base
    LF::lf.geometry
    LF::lf.fe
    LF::lf.io
    LF::lf.mesh.hybrid2d
    LF::lf.mesh
    LF::lf.mesh.utils
    LF::lf.uscalfe
)

add_executable(lecturecodes.convectiondiffusion.convergence convergence_main.cc cd_tools.cc ../../homeworks/UpwindQuadrature/mastersolution/upwindquadrature.cc)
target_compile_definitions(lecturecodes.convectiondiffusion.convergence PRIVATE CURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/${DIR}\")
target_compile_definitions(lecturecodes.convectiondiffusion.convergence PRIVATE CURRENT_BINARY_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
target_link_libraries(lecturecodes.convectiondiffusion.convergence
    PUBLIC Eigen3::Eigen
    LF::lf.assemble
    LF::lf.base
    LF::lf.geometry
    LF::lf.fe
    LF::lf.mesh.hybrid2d
    LF::lf.mesh
    LF::lf.mesh.utils
    LF::lf.uscalfe
    LF::lf.refinement
)