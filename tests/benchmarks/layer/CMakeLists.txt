if(NOT REAL)
    return()
endif(NOT REAL)
project(benchmark-layer)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-layer ${BIN})
set_tests_properties(benchmark-layer PROPERTIES LABELS slow)
