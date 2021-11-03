# define executable
set(sources
  eigensparseinit.cc
)

add_executable(lecturecodes.eigensparseinit ${sources})

target_link_libraries(lecturecodes.eigensparseinit
  PUBLIC Eigen3::Eigen
)
