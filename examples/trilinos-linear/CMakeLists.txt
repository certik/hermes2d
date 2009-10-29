if(NOT REAL)
    return()
endif(NOT REAL)
project(trilinos-linear)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
