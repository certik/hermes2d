enable_language(Fortran)
find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)

# examples
add_subdirectory(layer)
add_subdirectory(lshape)
add_subdirectory(smooth-iso)
add_subdirectory(smooth-aniso-x)
add_subdirectory(smooth-aniso-y)
add_subdirectory(singpert-aniso)
add_subdirectory(bessel)
add_subdirectory(screen)

