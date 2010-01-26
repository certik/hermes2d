if(NOT REAL)
    return()
endif(NOT REAL)
project(kellogg)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
