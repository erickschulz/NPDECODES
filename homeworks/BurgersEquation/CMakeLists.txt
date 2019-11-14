# Name this problem the same name as the containing folder.
get_filename_component(PROBLEM_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message(STATUS "Processing ${PROBLEM_NAME}")

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Provides a variable CURRENT_BINARY_DIR containing the path to the binary output directory.
add_definitions(-DCURRENT_BINARY_DIR=\"${CMAKE_BINARY_SOURCE_DIR}\")

# Provides build_problem(...) and build_test(...)
include(build_rules.cmake)

build_problem(${PROBLEM_NAME} mastersolution)
build_test(${PROBLEM_NAME} mastersolution)

build_problem(${PROBLEM_NAME} mysolution)
build_test(${PROBLEM_NAME} mysolution)
