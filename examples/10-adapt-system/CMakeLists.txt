if(NOT REAL)
    return()
endif(NOT REAL)
project(adapt-system)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
