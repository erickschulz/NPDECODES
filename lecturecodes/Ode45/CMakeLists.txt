# define executable
set(ode45sources ode45test.cc ode45.h)

add_executable(lecturecodes.ode45test ${ode45sources})

target_link_libraries(lecturecodes.ode45test
  PUBLIC Eigen3::Eigen
)
