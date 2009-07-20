if(NOT COMPLEX)
    return()
endif(NOT COMPLEX)
project(l2)

# use the complex version of the library:
set(HERMES ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
