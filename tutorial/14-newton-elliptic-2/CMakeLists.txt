if(NOT REAL)
    return()
endif(NOT REAL)
project(newton-elliptic-2)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
