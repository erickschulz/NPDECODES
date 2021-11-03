#!/usr/bin/env bash
# This script installs an up-2-date version of cmake on the travis build server.

# Exit immediately from this script upon error
set -e

# save current directory
dir=$(pwd)

mkdir -p ${DEPS_DIR} && cd ${DEPS_DIR}

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ] && [ ! -d "cmake-3.18.4" ]; then
  CMAKE_URL="https://github.com/Kitware/CMake/releases/download/v3.18.4/cmake-3.18.4-Linux-x86_64.tar.gz"
  mkdir cmake-3.18.4 && wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C cmake-3.18.4
elif [ "$(uname)" == "Darwin" ] && [ ! -d "cmake-3.18.4" ]; then
  CMAKE_URL="https://github.com/Kitware/CMake/releases/download/v3.18.4/cmake-3.18.4-Darwin-x86_64.tar.gz"
  mkdir cmake-3.18.4 && wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=3 -xz -C cmake-3.18.4
fi
export PATH=${DEPS_DIR}/cmake-3.18.4/bin:${PATH}

#Change back to where we left off.
cd $dir