cmake_minimum_required(VERSION 2.8.9)
project(ESKF)

IF (WIN32)
  set(CMAKE_CXX_FLAGS "-O2 -Wall")
ELSE()
  set(CMAKE_CXX_FLAGS "-Og -W -Wall -std=c++11")
ENDIF()

file(GLOB INCLUDES "include/*.h" "eigen/*.h" "eigen/src/*/*.h")
include_directories(./include eigen eigen/src/*/*)
file(GLOB SOURCES "src/*.cpp" "include/*.h" "eigen/*.h" "eigen/src/*/*.h")

add_executable(ESKF ${SOURCES})
