cmake_minimum_required(VERSION 3.16.4)
project(pso_uhpc LANGUAGES C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fopenmp")

add_subdirectory(src)
