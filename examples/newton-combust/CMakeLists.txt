if(NOT REAL)
    return()
endif(NOT REAL)
project(newton-combust)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
