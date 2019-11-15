# Provides the following targets in build folder (top-level only):

# test_solution_ON
# test_solution_OFF
# test_mastersolution
# test_mysolution
# templates

function(custom_targets HOMEWORKS DEVELOPERS LECTURECODES)

  set(DIR ${CMAKE_SOURCE_DIR})

  # Add custom targets for homeworks
  if(HOMEWORKS)
    add_custom_target(test_mastersolution
        COMMAND ${DIR}/scripts/run-tests.sh homeworks/*/*_test_mastersolution
    )
    add_custom_target(test_mysolution
        COMMAND ${DIR}/scripts/run-tests.sh homeworks/*/*_test_mysolution
    )
  endif()

  # Add custom targets for developers
  if(DEVELOPERS)
    add_custom_target(test_solution_ON
        COMMAND ${DIR}/scripts/run-tests.sh developers/*/*_test_solution_ON
    )
    add_custom_target(test_solution_OFF
        COMMAND ${DIR}/scripts/run-tests.sh developers/*/*_test_solution_OFF
    )
  endif()

  # Add custom target for creating templates from developers/ folder
  add_custom_target(templates
      COMMAND ${DIR}/scripts/make_templates.sh ${DIR}/developers/*/
  )

endfunction()
