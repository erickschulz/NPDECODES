# Provides variable PROBLEM_NAME
include(${CMAKE_SOURCE_DIR}/cmake/modules/build_variables.cmake)

# Provides functions build_problem and build_test
include(${CMAKE_SOURCE_DIR}/cmake/modules/build_rules.cmake)

set(MASTER_TARGET ${PROBLEM_NAME}_solution_${SOLUTION})
set(TEST_TARGET ${PROBLEM_NAME}_test_solution_${SOLUTION})

if(${SOLUTION})
  set(NAME_SUFFIX mastersolution)
else()
  set(NAME_SUFFIX mysolution)
endif()

message(STATUS "Processing ${PROBLEM_NAME}")
build_problem(${MASTER_TARGET} mastersolution ${PROBLEM_NAME}_${NAME_SUFFIX})
build_test(${TEST_TARGET} ${MASTER_TARGET} mastersolution ${PROBLEM_NAME}_test_${NAME_SUFFIX})
