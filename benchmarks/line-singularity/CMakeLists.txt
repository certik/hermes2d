if(NOT REAL)
    return()
endif(NOT REAL)
project(line-singularity)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
