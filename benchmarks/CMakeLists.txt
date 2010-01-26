enable_language(Fortran)
find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)

# examples
add_subdirectory(lshape)
add_subdirectory(layer)
add_subdirectory(smooth-iso)
add_subdirectory(smooth-aniso-x)
add_subdirectory(smooth-aniso-y)
add_subdirectory(singpert-aniso)
add_subdirectory(line-sing)
add_subdirectory(kellogg)

if(NOT WITH_TRILINOS)
  add_subdirectory(bessel)
  add_subdirectory(screen)
endif(NOT WITH_TRILINOS)

