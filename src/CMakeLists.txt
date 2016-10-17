cmake_minimum_required (VERSION 2.8.11)
project(TD1_PRCD)

add_library(cblas SHARED cblas.c)
add_library(perf SHARED perf.c)


if(DATA_TYPE_SINGLE)
  add_library(util util_float.c)
  add_executable(driver driver_single_precision.c)
  target_link_libraries(driver util cblas perf)
else(DATA_TYPE_SINGLE)
  add_library(util util.c)
  add_executable(driver driver_double_precision.c)
  target_link_libraries(driver util cblas perf)
endif(DATA_TYPE_SINGLE)