SET(CMAKE_CXX_FLAGS_ASAN "-O2 -g -fsanitize=address -fno-omit-frame-pointer")

SET(CMAKE_CXX_FLAGS_UBSAN "-O2 -g -fsanitize=undefined")
SET(CMAKE_EXE_LINKER_FLAGS_UBSAN "-fsanitize=undefined")
