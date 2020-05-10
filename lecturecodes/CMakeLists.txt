cmake_minimum_required(VERSION 3.10)

project(lecturecodes)

# hunter_add_package(Boost COMPONENTS math)
# find_package(Boost CONFIG REQUIRED math_c99)
# find_package(Boost CONFIG REQUIRED math_c99f)
# find_package(Boost CONFIG REQUIRED math_c99l)
# find_package(Boost CONFIG REQUIRED math_tr1)
# find_package(Boost CONFIG REQUIRED math_tr1f)
# find_package(Boost CONFIG REQUIRED math_tr1l)

## Subdirectories ##
# add_subdirectory(SimpleLinearFEM2D)  # needs MathGL
#add_subdirectory(EigenSparseMatrix)
# add_subdirectory(RatesOfConvergence) # needs MathGL
# Does not work until Boost math libraries are available
add_subdirectory(ConservativeFiniteVolumes)
