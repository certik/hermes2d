if(NOT REAL)
    return()
endif(NOT REAL)
project(trilinos-nonlinear)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
