IF(H2D_COMPLEX)

project(benchmark-bessel)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp bessel.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-bessel ${BIN})
set_tests_properties(benchmark-bessel PROPERTIES LABELS slow)

ENDIF(H2D_COMPLEX)
