# download dependencies
hunter_add_package(lehrfempp)
find_package(lehrfempp CONFIG REQUIRED)


# define executable
add_executable(hello main.cc)
# link executable with LehrFEM++ module lf::base
target_link_libraries(hello LF::lf.base) 