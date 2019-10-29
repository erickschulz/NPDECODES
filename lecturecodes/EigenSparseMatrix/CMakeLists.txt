# define executable
set(sources
  eigensparseinit.cc
)

add_executable(eigensparseinit ${sources})

target_link_libraries(eigensparseinit
  PUBLIC Eigen3::Eigen
)
