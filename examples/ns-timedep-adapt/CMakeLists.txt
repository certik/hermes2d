if(NOT REAL)
    return()
endif(NOT REAL)
project(ns-timedep-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
