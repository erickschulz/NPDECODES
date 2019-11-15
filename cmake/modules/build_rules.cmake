# Build rule for problems
function(build_problem TARGET DIR OUTPUT_NAME)
  # Defines SOURCES and LIBRARIES
  include(${DIR}/dependencies.cmake)

  add_executable(${TARGET} ${SOURCES})
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${OUTPUT_NAME})
  target_link_libraries(${TARGET}
	  ${LIBRARIES}
  )

  add_library(${TARGET}.shared SHARED ${SOURCES})
  set_target_properties(${TARGET}.shared PROPERTIES OUTPUT_NAME ${OUTPUT_NAME}.shared)
  target_link_libraries(${TARGET}.shared
	  ${LIBRARIES}
  )
endfunction(build_problem)

# Build rule for tests
function(build_test TARGET TARGET_TO_TEST DIR OUTPUT_NAME)
  # Defines SOURCES and LIBRARIES
  include(test/dependencies.cmake)
  include(GoogleTest)

  add_executable(${TARGET} ${SOURCES})
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${OUTPUT_NAME})
  target_link_libraries(${TARGET}
    ${LIBRARIES}
    ${TARGET_TO_TEST}.shared
  )

  gtest_discover_tests(${TARGET})
endfunction(build_test)
