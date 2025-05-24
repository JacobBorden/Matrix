cmake_minimum_required(VERSION 3.0.0)
project(bitmap VERSION 1.0.0)

include(CTest)
enable_testing()


add_library(matrix INTERFACE)
set_target_properties(matrix PROPERTIES LINKER_LANGUAGE CXX)
