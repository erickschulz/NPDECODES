# Call this script from the cmake build directory.
# It then goes through all subdirectories of ./developers/xxx and executes the command 
# - make -j2
# - ./xxx_test_mastersolution
# - make clean

set -e

for d in ./developers/*/ ;
do
  [[ $d =~ ./developers/(.*)/ ]]
  cmd="./${BASH_REMATCH[1]}_test_mastersolution";
  echo "-----------------------------------------------------------------------";
  echo "Compiling & Running tests in $d";
  echo "-----------------------------------------------------------------------";
  echo "Executing $cmd"
  (cd "$d" && make -j2 && eval $cmd && make clean);
done