# Provides variable PROBLEM_NAME
include(${CMAKE_SOURCE_DIR}/cmake/build_variables.cmake)

# Provides functions build_problem and build_test
include(${CMAKE_SOURCE_DIR}/cmake/build_rules.cmake)

function(build PROBLEM_NAME DIR)
  build_problem(${PROBLEM_NAME}_${DIR} ${DIR})
  build_test(${PROBLEM_NAME}_test_${DIR} ${PROBLEM_NAME}_${DIR} ${DIR})
endfunction(build)

message(STATUS "Processing ${PROBLEM_NAME}")
build(${PROBLEM_NAME} mastersolution)
build(${PROBLEM_NAME} mysolution)
