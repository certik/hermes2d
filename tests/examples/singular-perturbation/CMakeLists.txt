if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-singular-perturbation)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-singular-perturbation-projection ${BIN})
set_tests_properties(example-singular-perturbation-projection PROPERTIES LABELS slow)
add_test(example-singular-perturbation-coarse ${BIN} -coarse_mesh)
set_tests_properties(example-singular-perturbation-coarse PROPERTIES LABELS slow)
