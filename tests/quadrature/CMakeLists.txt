project(quadrature)

if(NOT UMFPACK_NO_BLAS)
	enable_language(Fortran)
	find_package(BLAS REQUIRED)
endif(NOT UMFPACK_NO_BLAS)
find_package(UMFPACK REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(quadrature-1 ${BIN})
