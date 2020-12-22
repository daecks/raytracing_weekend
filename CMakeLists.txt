cmake_minimum_required(VERSION 3.5)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -Wpedantic -O3)

# set the project name
project(RaytracingWeekend)

# add the executable
add_executable(RaytracingWeekend main.cpp)
