if(NOT REAL)
    return()
endif(NOT REAL)
project(lshape)

add_executable(${PROJECT_NAME} main.cpp)
include (../../examples/CMake.common)
