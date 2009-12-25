enable_language(Fortran)
find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)

# examples
add_subdirectory(layer)
add_subdirectory(lshape)
add_subdirectory(smooth)
add_subdirectory(anisopoly)
add_subdirectory(singpert)
add_subdirectory(bessel)
add_subdirectory(screen)

