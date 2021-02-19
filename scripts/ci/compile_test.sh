# Call this script from the cmake build directory.
# It then goes through all subdirectories of ./developers/xxx and executes the command 
# - make -j2
# - ./xxx_test_mastersolution
# - make clean

set -e

for d in ./developers/*/ ;
do
  [[ $d =~ ./developers/(.*)/ ]]
  cmd_test="./${BASH_REMATCH[1]}_test_mastersolution";
  cmd_mastersolution="./${BASH_REMATCH[1]}_mastersolution";
  
  if [[ $d =~ /CMakeFiles/ ]]; then
    continue
  fi
  echo "-----------------------------------------------------------------------";
  echo "Compiling & Running tests in $d";
  echo "-----------------------------------------------------------------------";
  cd "$d";
  make -j2;
  # Run unit tests if exists:
  if [[ -f "$cmd_test" ]]; then
    echo "Executing $cmd_test";
    eval $cmd_test
  else
    echo "*** WARNING: No unit tests found in $d ***";
  fi
  # Run solution if exists:
  if [[ -f "$cmd_mastersolution" ]]; then
    echo "Executing $cmd_mastersolution";
    if ! output=$(eval $cmd_mastersolution 2>&1) ; then
      echo "ERROR: ";
      printf "$output";
      exit 1;
    else
      echo "--> Success";
    fi
  else
    echo "*** WARNING: No mastersolution found in $d ***";
  fi
  
  make clean;
  cd ../..;
done
