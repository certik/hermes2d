if(NOT REAL)
    return()
endif(NOT REAL)
project(newton-timedep-flame)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
