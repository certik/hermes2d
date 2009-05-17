if(NOT REAL)
    return()
endif(NOT REAL)
project(multimesh)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
