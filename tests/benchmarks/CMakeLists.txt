find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
if(NOT UMFPACK_NO_BLAS)
	enable_language(Fortran)
	find_package(BLAS REQUIRED)
endif(NOT UMFPACK_NO_BLAS)

# examples
add_subdirectory(layer)
add_subdirectory(lshape)
add_subdirectory(smooth-iso)
add_subdirectory(smooth-aniso-x)
add_subdirectory(smooth-aniso-y)

if(NOT WITH_TRILINOS)
  add_subdirectory(bessel)
  add_subdirectory(screen)
endif(NOT WITH_TRILINOS)
