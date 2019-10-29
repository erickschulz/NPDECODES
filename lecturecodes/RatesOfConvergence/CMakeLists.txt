# define executable
set(sources
  estimated_order_convergence.h
  estimated_order_convergence.cc
)

add_executable(lecturecodes.rates_of_convergence  rates_of_convergence_main.cc ${sources})

target_link_libraries(lecturecodes.rates_of_convergence
  PUBLIC Eigen3::Eigen
  MathGL
)

# add_subdirectory(test)
