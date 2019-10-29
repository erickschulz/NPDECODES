cmake_minimum_required(VERSION 3.10)

## Hunter package manager ###
include("cmake/HunterGate.cmake")
HunterGate(
    URL "https://github.com/ruslo/hunter/archive/v0.23.104.tar.gz"
    SHA1 "92f60c3474c6d8749320f1b9ce4218a648493de3"
)

project(Lecturecodes)

# hunter_add_package(Boost COMPONENTS math)
# find_package(Boost CONFIG REQUIRED math_c99)
# find_package(Boost CONFIG REQUIRED math_c99f)
# find_package(Boost CONFIG REQUIRED math_c99l)
# find_package(Boost CONFIG REQUIRED math_tr1)
# find_package(Boost CONFIG REQUIRED math_tr1f)
# find_package(Boost CONFIG REQUIRED math_tr1l)

## Subdirectories ##
add_subdirectory(SimpleLinearFEM2D)
add_subdirectory(EigenSparseMatrix)
add_subdirectory(RatesOfConvergence)
# Does not work until Boost math libraries are available
# add_subdirectory(ConservativeFiniteVolumes)
