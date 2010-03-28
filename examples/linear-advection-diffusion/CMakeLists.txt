if(NOT REAL)
    return()
endif(NOT REAL)
project(linear-advection-diffusion)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
