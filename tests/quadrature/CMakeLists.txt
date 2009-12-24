project(quadrature)

enable_language(Fortran)
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(quadrature-1 ${BIN})
