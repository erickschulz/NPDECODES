# Provides variable PROBLEM_NAME
include(${CMAKE_SOURCE_DIR}/cmake/build_variables.cmake)

# Provides functions build_problem and build_test
include(${CMAKE_SOURCE_DIR}/cmake/build_rules.cmake)

set(MASTER_TARGET ${PROBLEM_NAME}_solution_${SOLUTION})
set(TEST_TARGET ${PROBLEM_NAME}_test_solution_${SOLUTION})

message(STATUS "Processing ${PROBLEM_NAME}")
build_problem(${MASTER_TARGET} mastersolution)
build_test(${TEST_TARGET} ${MASTER_TARGET} mastersolution)
