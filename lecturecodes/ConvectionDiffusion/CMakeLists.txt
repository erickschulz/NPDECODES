
add_custom_target(convectiondiffusion_run)


add_executable(lecturecodes.convectiondiffusion.layer layer_main.cc cd_tools.cc ../../homeworks/UpwindQuadrature/mastersolution/upwindquadrature.cc)
target_link_libraries(lecturecodes.convectiondiffusion.layer
    PUBLIC Eigen3::Eigen
    Boost::program_options
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


add_executable(lecturecodes.convectiondiffusion.convergence convergence.cc cd_tools.cc ../../homeworks/UpwindQuadrature/mastersolution/upwindquadrature.cc)
target_link_libraries(lecturecodes.convectiondiffusion.convergence
    PUBLIC Eigen3::Eigen
    Boost::program_options
    LF::lf.assemble
    LF::lf.base
    LF::lf.geometry
    LF::lf.fe
    LF::lf.io
    LF::lf.mesh.hybrid2d
    LF::lf.mesh
    LF::lf.mesh.utils
    LF::lf.uscalfe
    LF::lf.refinement
)


add_custom_target(lecturecodes.convectiondiffusion.layer_run
    COMMAND lecturecodes.convectiondiffusion.layer -o ${CMAKE_CURRENT_BINARY_DIR}/ && python3 ${CMAKE_CURRENT_SOURCE_DIR}/plot_layer.py ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/plot_layer.py
)
add_dependencies(convectiondiffusion_run lecturecodes.convectiondiffusion.layer_run)

add_custom_target(lecturecodes.convectiondiffusion.convergence_run
    COMMAND lecturecodes.convectiondiffusion.convergence -o ${CMAKE_CURRENT_BINARY_DIR}/ && python3 ${CMAKE_CURRENT_SOURCE_DIR}/plot_convergence.py ${CMAKE_CURRENT_BINARY_DIR}/results_errors.txt
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/plot_convergence.py
)
add_dependencies(convectiondiffusion_run lecturecodes.convectiondiffusion.convergence_run)
