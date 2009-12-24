if(NOT REAL)
    return()
endif(NOT REAL)
project(layer)

add_executable(${PROJECT_NAME} main.cpp)
include (../../examples/CMake.common)
