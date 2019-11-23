# Name this problem the same name as the containing folder.
get_filename_component(PROBLEM_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# Provides a variable CURRENT_SOURCE_DIR containing the path to this CMakeLists.txt file.
add_definitions(-DCURRENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Provides a variable CURRENT_BINARY_DIR containing the path to the binary output directory.
add_definitions(-DCURRENT_BINARY_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
