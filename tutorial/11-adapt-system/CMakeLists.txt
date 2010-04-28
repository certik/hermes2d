if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(adapt-system)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
