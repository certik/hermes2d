find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
if(NOT UMFPACK_NO_BLAS)
	enable_language(Fortran)
	find_package(BLAS REQUIRED)
endif(NOT UMFPACK_NO_BLAS)

# examples
add_subdirectory(newton-np-timedep-adapt-system)
add_subdirectory(iron-water)
add_subdirectory(saphir)
