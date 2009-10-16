if(NOT REAL)
    return()
endif(NOT REAL)
project(sea-breeze)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
