#! /bin/bash
# Creates mastersolution and templates from all-in-one C++ code
# Usage: make_templates <ProblemName>

script_dir=$(dirname $0)
root=$script_dir/..

for dir in $@; do
  ProblemName=$(basename $dir)
  echo "-- "$ProblemName

  input_dir=$root/developers/$ProblemName
  output_dir=$root/homeworks/$ProblemName

  # plain copy
  cp -r $input_dir/. $output_dir

  # problems
  cp -r $output_dir/mastersolution/. $output_dir/mysolution
  # Strip template parts 
  unifdef -DSOLUTION=1 -x 2 -m $output_dir/mastersolution/*.*
  # Strip solution 
  unifdef -DSOLUTION=0 -x 2 -m $output_dir/mysolution/*.*

  # tests
  if [ -d "$output_dir/mastersolution/test" ]; then
    unifdef -DSOLUTION=1 -x 2 -m $output_dir/mastersolution/test/*.*
    unifdef -DSOLUTION=0 -x 2 -m $output_dir/mysolution/test/*.*
  else
    echo "     Warning: Found no unit tests for $ProblemName"
  fi

  # templates it just plain copy of processed mysolution
  cp -r $output_dir/mysolution/. $output_dir/templates
done
