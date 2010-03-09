if(NOT REAL)
    return()
endif(NOT REAL)
project(thermoelasticity)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
