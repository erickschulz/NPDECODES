cmake_minimum_required(VERSION 3.10)

project(Developers)

# Set SOLUTION preprocessor directive.
option(SOLUTION "usage: cmake -DSOLUTION=ON/OFF .." ON)
if(SOLUTION)
  add_definitions(-DSOLUTION=true)
else()
  add_definitions(-DSOLUTION=false)
endif()
message(STATUS "Solution ${SOLUTION}")
# cmake -DCMAKE_CXX_FLAGS="-DSOLUTION=0" ..

add_subdirectory(BurgersEquation)
add_subdirectory(1DWaveAbsorbingBC)
