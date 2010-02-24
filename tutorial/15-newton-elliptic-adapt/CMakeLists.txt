if(NOT REAL)
    return()
endif(NOT REAL)
project(newton-elliptic-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
