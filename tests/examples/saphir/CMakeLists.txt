if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-saphir)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-saphir-projection ${BIN})
set_tests_properties(example-saphir-projection PROPERTIES LABELS slow)
add_test(example-saphir-coarse ${BIN} -coarse_mesh)
set_tests_properties(example-saphir-coarse PROPERTIES LABELS slow)
