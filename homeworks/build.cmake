# Provides variable PROBLEM_NAME
include(${CMAKE_SOURCE_DIR}/cmake/modules/build_variables.cmake)

# Provides functions build_problem and build_test
include(${CMAKE_SOURCE_DIR}/cmake/modules/build_rules.cmake)

function(build PROBLEM_NAME DIR)
  build_problem(${PROBLEM_NAME}_${DIR} ${DIR} ${PROBLEM_NAME}_${DIR})
  build_test(${PROBLEM_NAME}_test_${DIR} ${PROBLEM_NAME}_${DIR} ${DIR} ${PROBLEM_NAME}_test_${DIR})
endfunction(build)

message(STATUS "Processing ${PROBLEM_NAME}")
build(${PROBLEM_NAME} mastersolution)
build(${PROBLEM_NAME} mysolution)
