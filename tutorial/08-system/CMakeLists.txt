if(NOT REAL)
    return()
endif(NOT REAL)
project(system)

add_executable(${PROJECT_NAME} main.cpp)
include (../../examples/CMake.common)
