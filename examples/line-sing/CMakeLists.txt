if(NOT REAL)
    return()
endif(NOT REAL)
project(line-sing)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
