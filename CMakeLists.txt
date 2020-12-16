cmake_minimum_required(VERSION 3.5)

# Enable as many warnings as practical
#if(MSVC)
#    add_compile_options(/W4)
#else()
    add_compile_options(-Wall -Wextra -Wpedantic)
    #endif()

# set the project name
project(WeekendRayTracer)

# add the executable
add_executable(WeekendRayTracer main.cpp)
