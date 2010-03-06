if(NOT REAL)
    return()
endif(NOT REAL)
project(newton-timedep-heat-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
