# This script is called from .travis.yml
#
# It takes two inputs as environment variables:
#   COMPILER specifies the compiler to use (e.g. g++-7)
#   BUILD_TYPE specifies the CMAKE_BUILD_TYPE.
# 
# It builds NPDECODES using the provided compiler + cmake configuration.

# Exit immediately from this script upon error
set -e

# install new version of cmake:
mkdir -p ${HUNTER_ROOT}

# Install cmake
source $(dirname $0)/install_cmake.sh

# compile
cd ${TRAVIS_BUILD_DIR}
export CXX=${COMPILER}
$CXX --version

cmake -H. -BBuild -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DCMAKE_CXX_FLAGS=-g0 -Wdev
cd Build
make -j${NUM_PROC:-2}