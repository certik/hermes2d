if(NOT REAL)
    return()
endif(NOT REAL)
project(layer-2)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
