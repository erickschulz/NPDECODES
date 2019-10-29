cmake_minimum_required(VERSION 3.10)

#include("cmake/functions.cmake")
include("cmake/HunterGate.cmake")
include("cmake/sanitizers.cmake")
HunterGate(
    URL "https://github.com/ruslo/hunter/archive/v0.23.207.tar.gz"
    SHA1 "6e1b05a1f65c2aa492f7f43a43578e83ac62bbdd"
    LOCAL # use cmake/Hunter/config.cmake
)

project(NumPDE)

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

# Get Boost
hunter_add_package(Boost COMPONENTS chrono timer system filesystem program_options)
find_package(Boost CONFIG REQUIRED chrono timer system filesystem program_options)

# Problem: Actually boost::math is also needed, but there is a problem
# setting it up through hunter. 

# Get MathGL 
# Needs to be added to hunter first.
# hunter_add_package(MathGL)
# find_package(MathGL)

# Temporary MathGL fix
include(ExternalProject)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
message("-- Downloading MathGl")
set(MATHGL2_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/mathgl_install")
ExternalProject_Add(
    MathGL_external
    URL http://downloads.sourceforge.net/mathgl/mathgl-2.3.3.tar.gz
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/mathgl_binary
    DOWNLOAD_NO_PROGRESS 1
    CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
               "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
               -DOpenGL_GL_PREFERENCE=GLVND
               -DCMAKE_CXX_STANDARD=11
               -Denable-openmp=OFF
               -DMGL_HAVE_TYPEOF=0
               -DMGL_HAVE_C99_COMPLEX=0
               "-DMGL_LIB_INSTALL_DIR=${MATHGL2_INSTALL_DIR}/lib"
               "-DMGL_CGI_PATH=${MATHGL2_INSTALL_DIR}/share/mathgl"
               "-DCMAKE_INSTALL_PREFIX=${MATHGL2_INSTALL_DIR}"
		-DCMAKE_INSTALL_NAME_DIR=${MATHGL2_INSTALL_DIR}/lib
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/mathgl_install
)

set(MATHGL2_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/mathgl_install/include)
set(MATHGL2_LIBRARY_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}mgl${CMAKE_SHARED_LIBRARY_SUFFIX})
set(MATHGL2_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/mathgl_install/lib/${MATHGL2_LIBRARY_NAME}")
set(RPATH ${RPATH} "{CMAKE_CURRENT_BINARY_DIR}/mathgl_install/lib")

add_library(MathGL INTERFACE)
add_dependencies(MathGL MathGL_external)
target_include_directories(MathGL INTERFACE ${MATHGL2_INCLUDE_DIRS})
target_link_libraries(MathGL INTERFACE ${MATHGL2_LIBRARIES})

# Include Subdirectories
###############################################################################
add_subdirectory(homeworks)
add_subdirectory(lecturecodes)

add_custom_target(test_master
    COMMAND ${CMAKE_SOURCE_DIR}/scripts/run-tests.sh
)
