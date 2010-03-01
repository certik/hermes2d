if(NOT REAL)
    return()
endif(NOT REAL)
project(newton-timedep-ns)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
