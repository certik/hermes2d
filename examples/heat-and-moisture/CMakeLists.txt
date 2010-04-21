if(NOT REAL)
    return()
endif(NOT REAL)
project(heat-and-moisture)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
