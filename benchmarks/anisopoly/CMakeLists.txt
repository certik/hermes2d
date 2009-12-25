if(NOT REAL)
    return()
endif(NOT REAL)
project(anisopoly)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
