# define executable
set(ode45sources ode45test.cc)

add_executable(ode45test ${ode45sources})

target_link_libraries(ode45test
  PUBLIC Eigen3::Eigen
)
