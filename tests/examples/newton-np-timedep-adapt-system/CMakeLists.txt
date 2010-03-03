if(NOT REAL)
    return()
endif(NOT REAL)
project(example-np)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-np ${BIN})
