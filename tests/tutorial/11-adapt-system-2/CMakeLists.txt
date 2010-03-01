if(NOT REAL)
    return()
endif(NOT REAL)
project(tutorial-11-adapt-system-2)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-11-adapt-system-2 ${BIN})

