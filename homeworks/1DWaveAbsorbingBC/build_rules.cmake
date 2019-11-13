# Build rule for problems
function(build_problem PROBLEM_NAME DIR)
  # Defines SOURCES and LIBRARIES
  include(${DIR}/dependencies.cmake)

  add_executable(${PROBLEM_NAME}_${DIR} ${SOURCES})
  target_link_libraries(${PROBLEM_NAME}_${DIR}
	  ${LIBRARIES}
  )

  add_library(${PROBLEM_NAME}_${DIR}.shared ${SOURCES})
  target_link_libraries(${PROBLEM_NAME}_${DIR}.shared
	  ${LIBRARIES}
  )
endfunction(build_problem)

# Build rule for tests
function(build_test PROBLEM_NAME DIR)
  # Defines SOURCES and LIBRARIES
  include(test/dependencies.cmake)
  include(GoogleTest)

  add_executable(${PROBLEM_NAME}_test_${DIR} ${SOURCES})
  target_link_libraries(${PROBLEM_NAME}_test_${DIR}
    ${LIBRARIES}
  )

  gtest_discover_tests(${PROBLEM_NAME}_test_${DIR})
endfunction(build_test)
