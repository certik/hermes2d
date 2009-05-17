if(NOT REAL)
    return()
endif(NOT REAL)
project(newton-heat-timedep)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
