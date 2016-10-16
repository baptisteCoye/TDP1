cmake_minimum_required (VERSION 2.8.11)
project(TD1_PRCD)

add_library(util SHARED util.c)
add_library(cblas SHARED cblas.c)
add_library(perf SHARED perf.c)

add_executable(driver driver_double_precision.c)
target_link_libraries(driver util cblas perf)
