# This script installs an up-2-date version of cmake on the travis build server.

# Exit immediately from this script upon error
set -e

# save current directory
dir=$(pwd)

mkdir -p ${DEPS_DIR} && cd ${DEPS_DIR}

if [[ "${TRAVIS_OS_NAME}" == "linux" ]] && [ ! -d "cmake" ]; then
  CMAKE_URL="https://cmake.org/files/v3.11/cmake-3.11.1-Linux-x86_64.tar.gz"
  mkdir cmake && wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C cmake
elif [[ "${TRAVIS_OS_NAME}" == "osx" ]] && [ ! -d "cmake" ]; then
  CMAKE_URL="https://cmake.org/files/v3.11/cmake-3.11.1-Darwin-x86_64.tar.gz"
  mkdir cmake && wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=3 -xz -C cmake
fi
export PATH=${DEPS_DIR}/cmake/bin:${PATH}

#Change back to where we left off.
cd $dir