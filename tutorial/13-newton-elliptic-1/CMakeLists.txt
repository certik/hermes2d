if(NOT REAL)
    return()
endif(NOT REAL)
project(newton-elliptic-1)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
