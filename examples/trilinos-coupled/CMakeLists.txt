if(NOT REAL)
    return()
endif(NOT REAL)
project(trilinos-coupled)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
