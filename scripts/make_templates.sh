#! /bin/bash

script_dir=$(dirname $0)
root=$script_dir/..

for dir in $@; do

  ProblemName=$(basename $dir)
  echo "-- "$ProblemName

  input_dir=$root/developers/$ProblemName
  output_dir=$root/homeworks/$ProblemName

  cp -r $input_dir/. $output_dir
  cp -r $output_dir/mastersolution/. $output_dir/mysolution
  rename _mastersolution.cc _mysolution.cc $output_dir/test/*_mastersolution.cc
  cp $input_dir/test/*_mastersolution.cc $output_dir/test/

  unifdef -DSOLUTION=1 -x 2 -m $output_dir/mastersolution/* $output_dir/test/*_mastersolution.cc
  unifdef -DSOLUTION=0 -x 2 -m $output_dir/mysolution/* $output_dir/test/*_mysolution.cc

  cp -r $output_dir/mysolution/. $output_dir/templates

done
