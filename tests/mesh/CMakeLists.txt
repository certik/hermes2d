find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
if(NOT UMFPACK_NO_BLAS)
	enable_language(Fortran)
	find_package(BLAS REQUIRED)
endif(NOT UMFPACK_NO_BLAS)

# mesh tests
add_subdirectory(convert_to_quads)
add_subdirectory(convert_to_triangles)
add_subdirectory(refinements)
add_subdirectory(copy)
add_subdirectory(loader)

