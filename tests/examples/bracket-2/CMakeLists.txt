if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(bracket-2)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(bracket-2 ${BIN})

