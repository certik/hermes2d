IF(COMPLEX)

project(benchmark-bessel)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-bessel ${BIN})

ENDIF(COMPLEX)
