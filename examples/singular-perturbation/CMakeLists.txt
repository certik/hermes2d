if(NOT REAL)
    return()
endif(NOT REAL)
project(singpert-aniso)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
