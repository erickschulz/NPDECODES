#!/bin/bash

EXIT_CODE=0
FAILED_TESTS=()

for TEST_FILE in $@; do
  echo "Running $TEST_FILE"
  "$TEST_FILE"
  if [ $? -ne 0 ]; then
    PROBLEM_NAME=${TEST_FILE##*/}
    FAILED_TESTS+=("$PROBLEM_NAME")
    EXIT_CODE=1
  fi
done

#for TEST_FILE in $@; do
#    echo "Running $TEST_FILE"
#    "$TEST_FILE"
#    if [ $? -ne 0 ]; then
#        PROBLEM_NAME=${TEST_FILE##*/}
#        FAILED_TESTS+=("$PROBLEM_NAME")
#        EXIT_CODE=1
#    fi
#done

if [ $EXIT_CODE -ne 0 ]; then
  echo "The following tests have failed: "
  for PROBLEM_NAME in ${FAILED_TESTS[@]}; do
    echo " - $PROBLEM_NAME"
  done
  exit 1
else
  echo " ALL TESTS HAVE PASSED"
fi
