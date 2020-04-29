# Provides functions to create custom targets

# test targets
function(add_custom_test_target TARGET_NAME ARG)
  add_custom_target(${TARGET_NAME}
    COMMAND bash -c \"shopt -s nullglob &&  ${CMAKE_SOURCE_DIR}/scripts/run_tests.sh ${ARG}\"
  )
endfunction()

function(get_custom_test_targets)
  get_filename_component(DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  add_custom_test_target(test_${DIR}_mastersolution ${CMAKE_CURRENT_BINARY_DIR}/*/*_test_mastersolution)
  add_custom_test_target(test_${DIR}_mysolution ${CMAKE_CURRENT_BINARY_DIR}/*/*_test_mysolution)
endfunction()
