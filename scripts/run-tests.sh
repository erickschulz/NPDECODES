#!/bin/bash
EXIT_CODE=0
FAILED_TESTS=()
for e in homeworks/*/*test_mastersolution; do
  "$e"
  if [ $? -ne 0 ]; then
    FAILED_TESTS+=("$e")
    EXIT_CODE=1
  fi
done

for e in lecturecodes/*/test/*test; do
    "$e"
    if [ $? -ne 0 ]; then
        FAILED_TESTS+=("$e")
        EXIT_CODE=1
    fi
done


if [ $EXIT_CODE -ne 0 ]; then
  echo "The following tests have failed: "
  for e in $FAILED_TESTS; do
    echo " - $e"
  done
fi

exit $EXIT_CODE
