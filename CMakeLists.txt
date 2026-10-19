cmake_minimum_required(VERSION 3.16)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(projectName "SchedulabilityAnalysis")

project(${projectName})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_executable(${projectName} main.cpp)
target_compile_options(${projectName} PUBLIC -std=gnu++20)
