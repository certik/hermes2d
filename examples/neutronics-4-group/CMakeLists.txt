if(NOT REAL)
    return()
endif(NOT REAL)
project(neutronics-4-group)

set(HERMES ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
