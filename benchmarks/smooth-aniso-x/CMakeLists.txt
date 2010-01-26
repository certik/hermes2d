if(NOT REAL)
    return()
endif(NOT REAL)
project(smooth-aniso-x)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
