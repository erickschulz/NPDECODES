# Build rule for problems
function(build_problem TARGET DIR)
  # Defines SOURCES and LIBRARIES
  include(${DIR}/dependencies.cmake)

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
function(build_test TARGET TARGET_TO_TEST DIR)
  # Defines SOURCES and LIBRARIES
  include(test/dependencies.cmake)
  include(GoogleTest)

  add_executable(${TARGET} ${SOURCES})
  target_link_libraries(${TARGET}
    ${LIBRARIES}
    ${TARGET_TO_TEST}.shared
  )

  gtest_discover_tests(${TARGET})
endfunction(build_test)
