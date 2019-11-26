# Provides variable PROBLEM_NAME
include(${CMAKE_SOURCE_DIR}/cmake/modules/build_variables.cmake)

# Provides functions build_problem and build_test
include(${CMAKE_SOURCE_DIR}/cmake/modules/build_rules.cmake)

# pass correct arguemnts to build rules
function(build PROBLEM_NAME DIR SOLUTION)
  set(PROBLEM_TARGET ${PROBLEM_NAME}_${DIR})
  set(TEST_TARGET ${PROBLEM_NAME}_test_${DIR})

  # problem
  build_problem(${PROBLEM_TARGET}_dev mastersolution ${PROBLEM_TARGET})
  target_compile_definitions(${PROBLEM_TARGET}_dev PRIVATE SOLUTION=${SOLUTION})
  target_compile_definitions(${PROBLEM_TARGET}_dev.static PRIVATE SOLUTION=${SOLUTION})

  # tests
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mastersolution/test)
    build_test(${TEST_TARGET}_dev ${PROBLEM_TARGET}_dev mastersolution ${TEST_TARGET})
    target_compile_definitions(${TEST_TARGET}_dev PRIVATE SOLUTION=${SOLUTION})
  else()
    message(STATUS "*** Warning: Found no unit tests for ${PROBLEM_NAME} ***")
  endif()
endfunction(build)

# execute
message(STATUS "Processing ${PROBLEM_NAME}")

if(${MASTERSOLUTION})
  build(${PROBLEM_NAME} mastersolution 1)
endif()

if(${MYSOLUTION})
  build(${PROBLEM_NAME} mysolution 0)
endif()
