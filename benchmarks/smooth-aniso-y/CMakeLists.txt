if(NOT REAL)
    return()
endif(NOT REAL)
project(smooth-aniso-y)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
