cmake_minimum_required(VERSION 3.10)

include("cmake/HunterGate.cmake")
HunterGate(
    URL "https://github.com/cpp-pm/hunter/archive/v0.23.273.tar.gz"
    SHA1 "43ffd06cf0352babdbe28d333adb15379ab59681"
    LOCAL # use cmake/Hunter/config.cmake
)

project(NPDECODES)

message("Project source dir =" "${PROJECT_SOURCE_DIR}")
message("Project binary dir =" "${PROJECT_BINARY_DIR}")

# Check that this is not an in-source build:
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

#setup compiler
set(CMAKE_CXX_STANDARD 17)

#enable warnings
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
CHECK_CXX_COMPILER_FLAG(-Wextra COMPILER_SUPPORTS_WEXTRA)
CHECK_CXX_COMPILER_FLAG(-Wconversion COMPILER_SUPPORTS_WCONVERSION)

add_library(warning_flags INTERFACE)
if(COMPILER_SUPPORTS_WALL)
  target_compile_options(warning_flags INTERFACE -Wall)
endif()
if(COMPILER_SUPPORTS_WEXTRA)
  target_compile_options(warning_flags INTERFACE -Wextra)
endif()
if(COMPILER_SUPPORTS_WCONVERSION)
  target_compile_options(warning_flags INTERFACE -Wconversion)
endif()

# Get Dependencies
# (don't forget to update cmake/Config.cmake.in !
###############################################################################
# Get Eigen
hunter_add_package(Eigen)
find_package(Eigen3 CONFIG REQUIRED)

# Get LehrFEMpp
hunter_add_package(lehrfempp)
find_package(lehrfempp CONFIG REQUIRED)

# Get Google Test
hunter_add_package(GTest)
find_package(GTest CONFIG REQUIRED)
###############################################################################

# Add switches
option(HOMEWORKS "usage: cmake -DHOMEWORKS=ON/OFF .." ON)
message(STATUS "HOMEWORKS ${HOMEWORKS}")

option(DEVELOPERS "usage: cmake -DDEVELOPERS=ON/OFF .." ON)
message(STATUS "DEVELOPERS ${DEVELOPERS}")

option(LECTURECODES "usage: cmake -DLECTURECODES=ON/OFF .." ON)
message(STATUS "LECTURECODES ${LECTURECODES}")

option(MASTERSOLUTION "usage: cmake -DMASTERSOLUTION=ON/OFF .." ON)
message(STATUS "MASTERSOLUTION ${MASTERSOLUTION}")

option(MYSOLUTION "usage: cmake -DMYSOLUTION=ON/OFF .." ON)
message(STATUS "MYSOLUTION ${MYSOLUTION}")

# Add subdirectories
if(HOMEWORKS)
  add_subdirectory(homeworks)
endif()

if(DEVELOPERS)
  add_subdirectory(developers)
endif()

if(LECTURECODES)
  add_subdirectory(lecturecodes)
endif()
