# Build rule for problems
function(build_problem PROBLEM_NAME DIR)
  # Defines SOURCES and LIBRARIES
  include(${DIR}/dependencies.cmake)

  set(TARGET ${PROBLEM_NAME}_solution_${SOLUTION})

  add_executable(${TARGET} ${SOURCES})
  target_link_libraries(${TARGET}
	  ${LIBRARIES}
  )

  add_library(${TARGET}.shared ${SOURCES})
  target_link_libraries(${TARGET}.shared
	  ${LIBRARIES}
  )
endfunction(build_problem)

# Build rule for tests
function(build_test TARGET DIR)
  # Defines SOURCES and LIBRARIES
  include(test/dependencies.cmake)
  include(GoogleTest)

  set(TARGET ${PROBLEM_NAME}_test_solution_${SOLUTION})

  add_executable(${TARGET} ${SOURCES})
  target_link_libraries(${TARGET}
    ${LIBRARIES}
    ${PROBLEM_NAME}_solution_${SOLUTION}.shared
  )

  gtest_discover_tests(${TARGET})
endfunction(build_test)
