if(NOT REAL)
    return()
endif(NOT REAL)
project(example-bracket)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-bracket ${BIN})

