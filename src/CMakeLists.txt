cmake_minimum_required (VERSION 2.8.11)
project(TD1_PRCD)

add_library(util SHARED util.c)
add_library(ddot SHARED ddot.c)
add_library(perf SHARED perf.c)
add_library(dgemm SHARED dgemm.c)

add_executable(driver driver.c)
target_link_libraries(driver util ddot perf dgemm)
