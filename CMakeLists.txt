cmake_minimum_required(VERSION 2.8.11)
project(TD1_PRCD)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

option( DATA_TYPE_SINGLE "Set to ON if you want single precsion" OFF )


add_subdirectory(src)