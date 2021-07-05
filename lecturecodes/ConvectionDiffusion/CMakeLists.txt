

add_executable(lecturecodes.convectiondiffusion.main main.cc ../../homeworks/UpwindQuadrature/mastersolution/upwindquadrature.cc)

target_link_libraries(lecturecodes.convectiondiffusion.main
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


add_executable(lecturecodes.convectiondiffusion.convergence convergence.cc ../../homeworks/UpwindQuadrature/mastersolution/upwindquadrature.cc)

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
