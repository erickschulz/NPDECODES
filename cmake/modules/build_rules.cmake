# Build rule for problems
function(build_problem TARGET DIR OUTPUT_NAME)
  # Defines SOURCES and LIBRARIES
  include(${DIR}/dependencies.cmake)

  add_executable(${TARGET} ${SOURCES})
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${OUTPUT_NAME})
  target_compile_definitions(${TARGET} PRIVATE CURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/${DIR}\")
  target_compile_definitions(${TARGET} PRIVATE CURRENT_BINARY_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
  target_link_libraries(${TARGET}
	  ${LIBRARIES}
  )

  add_library(${TARGET}.static STATIC ${SOURCES})
  set_target_properties(${TARGET}.static PROPERTIES OUTPUT_NAME ${OUTPUT_NAME}.static)
  target_compile_definitions(${TARGET}.static PRIVATE CURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/${DIR}\")
  target_compile_definitions(${TARGET}.static PRIVATE CURRENT_BINARY_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
  target_link_libraries(${TARGET}.static
	  ${LIBRARIES}
  )
endfunction(build_problem)

# Build rule for tests
function(build_test TARGET TARGET_TO_TEST DIR OUTPUT_NAME)
  # Defines SOURCES and LIBRARIES
  include(${DIR}/test/dependencies.cmake)
  include(GoogleTest)

  add_executable(${TARGET} ${SOURCES})
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME ${OUTPUT_NAME})
  target_compile_definitions(${TARGET} PRIVATE CURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/test\")
  target_compile_definitions(${TARGET} PRIVATE CURRENT_BINARY_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
  target_link_libraries(${TARGET}
    ${LIBRARIES}
    ${TARGET_TO_TEST}.static
  )

  gtest_discover_tests(${TARGET})
endfunction(build_test)
