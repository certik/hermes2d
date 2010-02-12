find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
if(NOT UMFPACK_NO_BLAS)
	enable_language(Fortran)
	find_package(BLAS REQUIRED)
endif(NOT UMFPACK_NO_BLAS)

# adaptivity tests
add_subdirectory(cand_proj)
