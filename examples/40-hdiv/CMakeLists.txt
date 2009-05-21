if(NOT COMPLEX)
    return()
endif(NOT COMPLEX)
project(hdiv)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
