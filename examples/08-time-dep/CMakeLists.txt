if(NOT REAL)
    return()
endif(NOT REAL)
project(time-dep)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
