# Provides the following targets in build folder (top-level only):

# test_solution_ON
# test_solution_OFF
# test_mastersolution
# test_mysolution
# templates

function(custom_targets HOMEWORKS DEVELOPERS LECTURECODES)

  set(DIR ${CMAKE_SOURCE_DIR})

  # Add custom target for running all unit tests
  add_custom_target(test_mastersolution
      COMMAND shopt -s nullglob &&  ${DIR}/scripts/run-tests.sh homeworks/*/*_test_mastersolution developers/*/*_test_mastersolution
  )
  add_custom_target(test_mysolution
      COMMAND shopt -s nullglob && ${DIR}/scripts/run-tests.sh homeworks/*/*_test_mysolution developers/*/*_test_mysolution
  )

  # Add custom target for creating templates from developers/ folder
  add_custom_target(templates
      COMMAND ${DIR}/scripts/make_templates.sh ${DIR}/developers/*/
  )

endfunction()
