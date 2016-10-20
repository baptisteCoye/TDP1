cmake_minimum_required(VERSION 2.8.11)
project(TD1_PRCD)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

option( DATA_TYPE_SINGLE "Set to ON if you want single precsion" OFF )


add_subdirectory(src)
