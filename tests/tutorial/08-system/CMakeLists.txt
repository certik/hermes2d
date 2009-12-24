if(NOT REAL)
    return()
endif(NOT REAL)
project(tutorial-08-system)

add_executable(${PROJECT_NAME} main.cpp)
include (../../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-08-system ${BIN})

