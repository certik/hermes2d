if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(benchmark-layer)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-layer-projection ${BIN})
set_tests_properties(benchmark-layer-projection PROPERTIES LABELS slow)
add_test(benchmark-layer-coarse ${BIN} -coarse_mesh)
set_tests_properties(benchmark-layer-coarse PROPERTIES LABELS slow)
