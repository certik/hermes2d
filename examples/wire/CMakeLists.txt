if(NOT COMPLEX)
    return()
endif(NOT COMPLEX)


project(wire)

set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
